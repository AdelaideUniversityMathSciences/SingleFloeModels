function [ colour ] = ColourCode( range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Colour code
Colour = linspace(0,1,length(range));
colour = cell(length(Colour),1);
for n = 1:length(Colour)
 C = Colour(n);
 colour(n) = {[0 1-C C]};
end

end


