function w = calculate_width(data,diff)

a = data > max (data)/2;
first = find(a, 1, 'first');
final = find(a, 1, 'last');

w = diff*(final-first);

end