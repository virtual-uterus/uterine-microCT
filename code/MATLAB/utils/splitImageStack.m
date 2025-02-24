function stack_cell = splitImageStack(img_stack, split_nb)
%SPLITIMAGESTACK Splits an image stack into three parts (top, left, and
%right) to be processed indiviudally.
%   
%   Input:
%    - img_stack, stack of images to split.
%    - split_nb,  slice at which to split into left and right.

%   Return:
%    - stack_cell, cell array containing top, left, and right stacks
%       (in that order) 
if split_nb < 1
    error("The split number is smaller than the number of slices in the stack");

elseif split_nb > size(img_stack, 3)
    error("The split number is greater than the number of slices in the stack");
end

top_stack = img_stack(:, :, 1:split_nb-1);
left_stack = img_stack(:, 1:round(end/2), split_nb:end);
right_stack = img_stack(:, round(end/2)+1:end, split_nb:end);

stack_cell = {top_stack, left_stack, right_stack};