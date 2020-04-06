struct UnionJack <: Method
    axis_method::SOS2Method
end

# Logarithmic is default method for x_1 and x_2 axis sos2 constraints.
UnionJack() = UnionJack(Logarithmic())
