def box_muller(num, mean, stddev)
  samples = []
  loop do
    theta = 2 * Math::PI * Kernel.rand
    rho = Math.sqrt(-2 * Math.log(1 - Kernel.rand))
    x = mean + stddev * rho
    samples << x * Math.cos(theta) if samples.length < num
    samples << x * Math.sin(theta) if samples.length < num
    return samples if samples.length >= num
  end
end
