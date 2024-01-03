function dark_channel = get_dark_channel(img, win_size)

[m, n, ~] = size(img);



pad_size = floor(win_size/2);

padded_image = padarray(img, [pad_size pad_size], Inf);

dark_channel = zeros(m, n); 

for j = 1 : m
    for i = 1 : n
        patch = padded_image(j : j + (win_size-1), i : i + (win_size-1), :);
        dark_channel(j,i) = min(patch(:));
     end
end

end