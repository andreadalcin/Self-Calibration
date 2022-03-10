function [mrkrSize, syms, colorScheme, colorOutlier, imageAlpha, mrkrFaceAlpha] = parseOptsClustDisplay(opts)
%PARSEOPTSDISPLAY


% markerSize
if(isfield(opts,'mrkrSize'))
    mrkrSize = opts.mrkrSize;
else
    mrkrSize = 30;
end

% markerFaceAlpha
if(isfield(opts,'mrkrFaceAlpha'))
    mrkrFaceAlpha = opts.mrkrFaceAlpha;
else
    mrkrFaceAlpha = 1;
end

% symbols
if(isfield(opts,'syms'))
    syms = opts.syms;
else
    syms = 'osd<v>^ph';
end

% color scheme
if(isfield(opts,'scheme'))
    colorScheme = opts.scheme;
else
    colorScheme = 'Accent';
end

%color outlier
if(isfield(opts,'colorOutlier'))
    colorOutlier = opts.colorOutlier;
else
    colorOutlier = [0.3,0.3,0.3];
end

%alpha 
if(isfield(opts,'imageAlpha'))
    imageAlpha = opts.imageAlpha;
else
    imageAlpha = 0.9;
end




end

