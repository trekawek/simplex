require 'matrix'
require 'mathn'
require 'rational'

class Vector
  def < value
    to_a.each { |v| return false if v >= value }
    true
  end
end

class Matrix
  # Some ugly function that generates vector with values inserted
  # into indices. If ident != nil then on 1 is set on the ident
  # place.
  def self.prepare_from_array(values, indices, size, ident = nil)
    a = Array.new(size, 0)
    values.each_index { |i| a[indices[i]] = values[i] }
    a[((0..(size-1)).to_a - indices).sort[ident]] = 1 if ident
    Matrix[a].transpose
  end

  def < value
    row_vectors.each do |row|
      row.to_a.each { |v| return false if v >= value }
      true
    end
  end

  def get_columns(*indices)
    column_vectors.values_at(*indices)
  end

  def get_matrix_from_columns(*indices)
    Matrix.columns(get_columns(*indices))
  end

  def get_matrix_from_column(col, *indices)
    Matrix[column(col).to_a.values_at(*indices)]
  end
end

class Simplex
  def initialize n_M, b_M, c_V, b_V
    @n_M = n_M
    @b_M = b_M
    @c_V = c_V
    @b_V = b_V
  end

  def result
    return @result if @result
    @result = simplex @n_M, @b_M, @c_V, @b_V, false
  end

  def print_solution
    @result = simplex @n_M, @b_M, @c_V, @b_V, true
    puts
    puts "Rezultat: #{@result[:result]}"
    puts "Wynik: #{@result[:point].inspect}" if @result[:result] == :ok
    puts "Wartosc funkcji: #{get_value}"
    true
  end

  def get_value
    r = result
    return nil unless r[:result] == :ok
    v = 0
    c = @c_V.column(0).to_a
    r[:point].each_index { |i| v += r[:point][i] * c[i] }
    v
  end

  private 
  def create_matrix_A n_M, b_M
    Matrix.columns(n_M.column_vectors + b_M.column_vectors)
  end

  def create_first_x b_V, b_indices, a_M
    Matrix.prepare_from_array b_V.column(0).to_a, b_indices, a_M.column_size
  end

  def find_j alpha
    alpha_array = alpha.row(0).to_a
    alpha_array.index(alpha_array.max)
  end

  def simplex n_M, b_M, c_V, b_V, debug = false
    return @result if @result and !debug
    a_M = create_matrix_A n_M, b_M
    b_indices = ((n_M.column_size)..(a_M.column_size - 1)).to_a
    x_ = create_first_x b_V, b_indices, a_M
    history = [x_.column(0).to_a.join(' ')]
    while true
      n_indices = (0..(a_M.column_size - 1)).to_a - b_indices
      b_M = a_M.get_matrix_from_columns(*b_indices)
      n_M = a_M.get_matrix_from_columns(*n_indices)
      b_inversed_M = b_M.inverse
      if debug
        puts "B: " + b_M.inspect
        puts "B^-1: " + b_inversed_M.inspect
        puts "N: " + n_M.inspect
      end

      c_B = c_V.get_matrix_from_column(0, *b_indices)
      c_N = c_V.get_matrix_from_column(0, *n_indices)
      alpha = c_B * b_inversed_M * n_M - c_N
      puts "Alpha: " + alpha.inspect if debug
      return :result => :ok, :point => x_.column(0).to_a if alpha < 0
      
      j = find_j(alpha)
      a_j = a_M.column(j)
      y_j = b_inversed_M * a_j
      puts "y_#{j}: " + y_j.inspect if debug
      return :result => :no_solution if y_j < 0

      b_ = (b_inversed_M * b_V).column(0).to_a
      l = []
      y_j.to_a.each_index { |i| l << b_[i] / y_j[i] if y_j[i] > 0 }
      l = l.min
      return :result => :ok, :point => Array.new(x_.row_size, 0) if l == 0

      v_j = Matrix.prepare_from_array(
        (b_inversed_M * a_j * -1).to_a,
        b_indices,
        a_M.column_size,
        j)
      puts "v_#{j}: " + v_j.inspect if debug

      x_ = x_ + l * v_j
      puts "new x: " + x_.inspect if debug

      x_a = x_.column(0).to_a
      x_s = x_a.join(' ')
      return :result => :loop if history.index(x_s)
      history << x_s

      b_indices = []
      x_a.each_index { |i| b_indices << i if x_a[i] != 0 }

      puts
    end
  end
end

n_M = Matrix[[4, 4], [1, 0], [0, 1]]
b_M = Matrix.identity(3)
c_V = Matrix[[-1, -2, 0, 0, 0]].transpose
b_V = Matrix[[12, 2, 2]].transpose

s = Simplex.new n_M, b_M, c_V, b_V
s.print_solution
