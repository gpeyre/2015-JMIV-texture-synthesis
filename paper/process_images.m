outdir = './used-images/';
if not(exist(outdir))
    mkdir(outdir);
end

a = 'imageslist.txt';
fid = fopen(a, 'rt');
tline = fgetl(fid);
while ischar(tline)
    i = findstr(tline, '/');
    u = [outdir tline(1:i(end))];
    if not(exist(u))
        mkdir(u);
    end
    copyfile(tline, [outdir tline], 'f');
    tline = fgetl(fid);
end