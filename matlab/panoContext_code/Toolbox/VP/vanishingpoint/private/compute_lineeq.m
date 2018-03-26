function lines = compute_lineeq(lines)

for i=1:length(lines)
	
	lineeq = line_equation_from_two_points(lines(i).point1, lines(i).point2);
	
	if lines(i).lineclass == 1
		if distance_of_point_to_line(lineeq, [0 -10000000]) < 0
			lineeq = -lineeq;
		end
	else
		if distance_of_point_to_line(lineeq, [10000000 0]) < 0
			lineeq = -lineeq;
		end
	end
	
	lines(i).lineeq = lineeq;
end
