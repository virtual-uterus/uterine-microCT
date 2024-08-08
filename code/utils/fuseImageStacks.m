function img_stack = fuseImageStacks(stack_cell)
%FUSEIMAGESTACKS Fuses the image stacks from a stack cell into a single
%stack of images.
%   
%   Input:
%    - stacks_cell, cell array containing top, left, and right stacks
%       (in that order) 
%
%   Return:
%    - img_stack, stack of fused images.
top_stack = stack_cell{1};
left_stack = stack_cell{2};
right_stack = stack_cell{3};

if ~isequal(size(left_stack, 1), size(right_stack, 1)) || ~isequal(size(left_stack, 3), size(right_stack, 3))
    error("The dimensions of the left and right stacks do not agree")

elseif size(top_stack, 2) ~= size(right_stack, 2) + size(left_stack, 2)
    error("The width of the three stacks do not agree")

end

img_stack = cat(3, top_stack, cat(2, left_stack, right_stack));