require 'matrix'
require 'mathn'
require 'rational'
require 'rubygems'
require 'terminal-table/import'

class Vector
  def >= v
    v = v.to_a
    to_a.each_index do |i|
      return true if self[i] > v[i]
      return false if self[i] < v[i]
    end
    true
  end

  def > v
    self >= v and self != v
  end

  def < v
    not self >= v
  end

  def <= v
    self < v or self == v
  end
end

class Gomory
  def initialize constants, tab
    @table = GomoryTable.new constants, tab
  end

  def print_solution
    i = 1
    tab = @table
    while true
      puts "\n ###### Iteration #{i} ###### \n"
      puts "Table:\n#{tab}"
      j, k = tab.get_jk
      break if !j
      puts "j = #{j}, k = #{k}"
      e = tab.get_e j, k
      puts "e = #{e.inspect}"
      l = tab.get_lambda j, e
      puts "lambda = #{l}"
      tab.add_new_x j, l
      puts "With new variable:\n#{tab}"
      tab = tab.sub_with_last k
      puts "After substitute:\n#{tab}"
      i += 1
    end
  end

  def result
    tab = @table
    while true
      j, k = tab.get_jk
      break if !j
      e = tab.get_e j, k
      l = tab.get_lambda j, e
      tab.add_new_x j, l
      tab = tab.sub_with_last k
    end
    tab
  end
end

class GomoryTable
  def initialize constants, tab, col_variables = nil, row_variables = nil
    @constants = constants
    @table_matrix = Matrix[*tab]
    cols = tab[0].length
    rows = constants.length
    @row_variables ||= ['x0'] + (((cols + 1)..(cols + rows - 1)).to_a.map { |i| "x#{i}" })
    @col_variables ||= (1..(cols)).map { |i| "-x#{i}" }
  end

  def get_jk
    min = @constants.min
    return nil if min >= 0
    j = @constants.index(min)
    max_v = nil; k = nil
    columns = @table_matrix.column_vectors
    columns.each_index do |i|
      v = columns[i]
      next if v[j] >= 0
      if !max_v or max_v < v
        max_v = v
        k = i
      end
    end
    return j, k
  end

  def get_e j, k
    e = []
    columns = @table_matrix.column_vectors
    0.upto(@table_matrix.column_size - 1) do |i|
      e[i] = 1 and next if i == k
      v1 = columns[i] * -1
      v2 = columns[k] * -1
      f = 1
      while true
        f += 1 and next if v1 * (1 / f) >= v2
        e[i] = f - 1 and break
      end
    end
    e
  end

  def get_lambda j, e
    row = @table_matrix.row(j).to_a
    set = []
    row.each_index { |i| set << -row[i] / e[i] if row[i] < 0 }
    set.max
  end

  def add_new_x j, l
    row = @table_matrix.row(j).to_a
    c = (@constants[j] / l).floor
    new_row = []
    row.each_index { |i| new_row << (row[i] / l).floor }
    @constants << c
    @row_variables << "x#{@constants.length + 1}"
    @table_matrix = Matrix[*(@table_matrix.row_vectors.map { |v| v.to_a } + [new_row])]
  end

  def sub_with_last k
    c = @constants.last
    new_const = []
    new_rows = []
    row = @table_matrix.row(@constants.length - 1).to_a
    row.each_index do |i|
      if i == k
        row[i] = -1
      else
        row[i] *= -1
      end
    end

    rows = @table_matrix.row_vectors
    new_rows = []
    new_constants = []
    0.upto(rows.length - 2) do |i|
      m = rows[i][k] * -1
      new_c = @constants[i] - c * m
      new_row = []
      rows[i].to_a.each_index do |j|
        if j == k
          new_row[j] = row[j] * m
        else
          new_row[j] = rows[i][j] + row[j] * m
        end
      end
      new_rows << new_row
      new_const << new_c
    end
    GomoryTable.new new_const, new_rows
  end

  def to_s
    my_table = table do |t|
      t.headings = ['', 'const'] + @col_variables
      rows = @table_matrix.row_vectors
      rows.each_index { |i| t << [@row_variables[i], @constants[i]] + rows[i].to_a }
    end
    my_table.to_s
  end
end

g = Gomory.new [0, -8, -6], [[-2, -3], [-1, -4], [-3, -1]]
g.print_solution
