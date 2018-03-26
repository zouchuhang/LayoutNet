function tex_responses = compute_texton_response(im, tim)

if isempty(tim)
    tex_responses = [];
    return;
end

tex_responses = zeros(size(im, 1), size(im, 2), numel(tim));

fims = fbRun(tim, im);

for i = 1:length(fims)
    tex_responses(:, :, i) = abs(fims{i});
end
