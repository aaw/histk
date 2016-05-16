class LosslessHistogram
  def initialize(values=[])
    @values = values.sort
  end

  def quantile(q)
    index = (q * @values.length).round
    index = 0 if index < 0
    index = -1 if index >= @values.length
    @values[index]
  end

  def count(v)
    @values.find_index{ |x| x > v } || @values.length
  end
end
