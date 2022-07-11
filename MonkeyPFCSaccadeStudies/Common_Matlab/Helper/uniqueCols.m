function y = uniqueCols(x)

[srtX, srtIdx] = sortrows(x');
dX = diff(srtX, 1, 1);
unqIdx = [true; any(dX, 2)];
y = x(:,srtIdx(unqIdx));