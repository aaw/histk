require 'redis'
require 'test/unit'

require_relative 'lossless_histogram'
require_relative 'box_muller'

def restart_redis(args=nil)
  `/bin/bash /opt/histk/restart_redis.sh #{args}`
end

def rm_redis_file(conn, filename)
  rdb_file = [conn.config(:get, 'dir')['dir'], filename].join('/')
  `rm -f #{rdb_file}`
end

class TestHistk < Test::Unit::TestCase
  def setup
    @conn = Redis.new(host: 'localhost', port: ENV['REDIS_PORT'])
    @r = @conn.client
    rm_redis_file(@conn, 'dump.rdb')
    rm_redis_file(@conn, 'appendonly.aof')
    restart_redis
  end

  def test_empty
    exception = assert_raise(Redis::CommandError) do
      @r.call(%w(histk.quantile s 0.5))
    end
    err = 'ERR empty histogram.'
    assert_equal(err, exception.message)
  end

  def test_type
    @r.call(%w(histk.add s 1))
    assert_equal('aaw-histk', @r.call(%w(type s)))
  end

  def test_wrong_type_on_histk_key
    @r.call(%w(histk.add a 1))
    exception = assert_raise(Redis::CommandError) do
      @conn.incr('a')
    end
    err = 'WRONGTYPE Operation against a key holding the wrong kind of value'
    assert_equal(err, exception.message)
  end

  def test_wrong_type_on_non_histk_key
    @conn.set('a', 100)
    [%w(histk.add a 1), %w(histk.count a), %w(histk.quantile a 0.5),
     %w(histk.mergestore b a), %w(histk.mergestore a c),
     %w(histk.resize a 8)].each do |cmd|
      exception = assert_raise(Redis::CommandError) do
        @r.call(cmd)
      end
      err = 'WRONGTYPE Operation against a key holding the wrong kind of value'
      assert_equal(err, exception.message)
    end
  end

  def test_basic_resize
    [16,2,8,32,128,1024].each do |size|
      assert_equal(size, @r.call(['histk.resize', 's', size]))
    end
  end

  def test_huge_resize
    exception = assert_raise(Redis::CommandError) do
      @r.call(%w(histk.resize s 2049))
    end
    err = 'ERR invalid size: number of centroids must be at most 2048.'
    assert_equal(err, exception.message)
  end

  def test_add
    assert_equal(1, @r.call(%w(histk.add s 100)))
    assert_equal(2, @r.call(%w(histk.add s 100)))
    assert_equal(3, @r.call(%w(histk.add s 200)))
    assert_equal(4, @r.call(%w(histk.add s 200)))
    assert_equal(5, @r.call(%w(histk.add s 300)))
    assert_equal(6, @r.call(%w(histk.add s 400)))
    assert_equal(6, @r.call(%w(histk.count s)))
    assert_equal('200', @r.call(%w(histk.quantile s 0.5)))
  end

  def test_add_multi
    assert_equal(6, @r.call(['histk.add', 's', 100, 2, 200, 2, 300, 1, 400, 1]))
    assert_equal(6, @r.call(%w(histk.count s)))
    assert_equal('200', @r.call(%w(histk.quantile s 0.5)))
  end

  def test_count
    @r.call(%w(histk.resize s 8))
    count = 0
    (1..10).each do |i|
      (1..50).each do |j|
        @r.call(['histk.add', 's', j, i])
        count += i
        assert_equal(count, @r.call(%w(histk.count s)))
      end
    end
  end

  def test_quantile_too_large
    exception = assert_raise(Redis::CommandError) do
      @r.call(%w(histk.quantile s 1.1))
    end
    err = 'ERR argument must be in the range [0.0, 1.0].'
    assert_equal(err, exception.message)
  end

  def test_quantile_too_small
    exception = assert_raise(Redis::CommandError) do
      @r.call(%w(histk.quantile s -0.1))
    end
    err = 'ERR argument must be in the range [0.0, 1.0].'
    assert_equal(err, exception.message)
  end

  def test_quantile_uniform
    @r.call(%w(histk.resize s 8))
    (1..100).each do |i|
      @r.call(['histk.add', 's', i])
    end
    assert_equal('1', @r.call(%w(histk.quantile s 0.0)))
    assert_equal('100', @r.call(%w(histk.quantile s 1.0)))
    abs_error = 1.5  # Arbitrary.
    (1..9).each do |i|
      actual = @r.call(['histk.quantile', 's', (0.1 * i).round(1)])
      assert_operator((10 * i - actual.to_f).abs, :<, abs_error)
    end
  end

  def test_count_uniform
    @r.call(%w(histk.resize s 32))
    (1..100).each do |i|
      @r.call(['histk.add', 's', i])
    end
    assert_equal(0, @r.call(%w(histk.count s 0.0)))
    assert_equal(100, @r.call(%w(histk.count s 100.0)))
    abs_error = 1.5  # Arbitrary.
    (1..100).each do |i|
      actual = @r.call(['histk.count', 's', i])
      assert_operator((i - actual.to_f).abs, :<, abs_error)
    end
  end

  def test_quantile_normal
    Kernel.srand(1234)  # Arbitrary fixed seed for reproducibility.
    values = box_muller(10000, 0, 1)
    values.each do |sample|
      @r.call(['histk.add', 's', sample])
    end
    h = LosslessHistogram.new(values)
    expected_error = 0.1  # Arbitrary
    [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999].each do |q|
      expected = h.quantile(q)
      actual = @r.call(['histk.quantile', 's', q]).to_f
      error = (expected - actual).abs
      assert_operator((expected - actual).abs, :<, expected_error)
    end
  end

  def test_count_normal
    Kernel.srand(1234)  # Arbitrary fixed seed for reproducibility.
    values = box_muller(10000, 0, 1)
    values.each do |sample|
      @r.call(['histk.add', 's', sample])
    end
    h = LosslessHistogram.new(values)
    expected_error = 0.02  # Arbitrary
    [-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4].each do |v|
      expected = h.count(v)
      actual = @r.call(['histk.count', 's', v]).to_f
      next if expected == 0 && actual == 0
      error = (expected - actual).abs / [expected.abs, actual.abs].max
      assert_operator(error, :<, expected_error)
    end
  end

  def test_ping_dataset
    values = []
    File.open('/opt/histk/ping_data').each_line do |line|
      @r.call(['histk.add', 's', line.strip])
      values << line.strip.to_f
    end
    h = LosslessHistogram.new(values)
    low = [0.00001, 0.0001, 0.001, 0.01, 0.1]
    mid = [0.25, 0.5, 0.75]
    high = [0.9, 0.99, 0.999, 0.9999, 0.99999]
    (low + mid + high).each do |q|
      expected = h.quantile(q)
      actual = @r.call(['histk.quantile', 's', q]).to_f
      relative_error = (expected - actual).abs / [expected.abs, actual.abs].max
      assert_operator(relative_error, :<, 0.05)
    end
  end

  def test_merge
    (1..100).each { |i| @r.call(['histk.add', 's', i]) }
    (101..200).each { |i| @r.call(['histk.add', 't', i]) }
    (201..300).each { |i| @r.call(['histk.add', 'u', i]) }
    assert_equal(200, @r.call(%w(histk.mergestore v s t)))
    assert_equal(300, @r.call(%w(histk.mergestore w s t u)))
    error = 1.5  # Arbitrary
    x = @r.call(%w(histk.quantile s 0.50)).to_f
    y = @r.call(%w(histk.quantile v 0.25)).to_f
    z = @r.call(%w(histk.quantile w 0.17)).to_f
    assert_operator((x - y).abs, :<, error)
    assert_operator((y - z).abs, :<, error)
  end

  def test_rdb
    @r.call(%w(histk.add s 100))
    @r.call(%w(histk.add s 200))
    @r.call(%w(histk.add s 300))
    args = [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    @r.call(['save'])
    restart_redis
    new_qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    assert_equal(qs, new_qs)
  end

  def test_aof
    restart_redis '--appendonly yes --appendfsync always'
    @r.call(%w(histk.add s 100))
    @r.call(%w(histk.add s 200))
    @r.call(%w(histk.add s 300))
    args = [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    restart_redis '--appendonly yes'
    new_qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    assert_equal(qs, new_qs)
  end

  def test_rewrite_aof
    restart_redis '--appendonly yes'
    @r.call(%w(histk.add s 100))
    @r.call(%w(histk.add s 100))
    @r.call(%w(histk.add s 200))
    @r.call(%w(histk.add s 100))
    10.times { @r.call(%w(histk.add s 200)) }
    @r.call(%w(histk.add s 300))
    args = [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    @r.call(['bgrewriteaof'])
    while @conn.info('persistence')['aof_rewrite_scheduled'] != '0' && \
          @conn.info('persistence')['aof_rewrite_in_progress'] != '0'
      sleep 0.1
    end
    aof_rewrite_status = @conn.info('persistence')['aof_last_bgrewrite_status']
    assert_not_equal(aof_rewrite_status, 'err')
    restart_redis '--appendonly yes'
    new_qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    assert_equal(qs, new_qs)
  end

  def test_rewrite_aof_with_resize
    restart_redis '--appendonly yes'
    @r.call(%w(histk.resize s 4))
    (1..100).each { |i| @r.call(['histk.add', 's', i]) }
    args = [0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    qs = args.map{ |q| @r.call(['histk.quantile', 's', q]) }
    @r.call(%w(histk.resize t 4))
    (1..50).each { |i| @r.call(['histk.add', 't', i]) }
    @r.call(['bgrewriteaof'])
    while @conn.info('persistence')['aof_rewrite_scheduled'] != '0' && \
          @conn.info('persistence')['aof_rewrite_in_progress'] != '0'
      sleep 0.1
    end
    aof_rewrite_status = @conn.info('persistence')['aof_last_bgrewrite_status']
    assert_not_equal(aof_rewrite_status, 'err')
    restart_redis '--appendonly yes'
    # t is empty but should have retained its size (4 centroids)
    (51..100).each { |i| @r.call(['histk.add', 't', i]) }
    new_qs = args.map{ |q| @r.call(['histk.quantile', 't', q]) }
    assert_equal(qs, new_qs)
  end
end
