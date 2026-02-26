function [spec] = fSpecImport(filename,d,scaleUp,sumA)

% Pre-porcesses the input wave spectrum

% filename - path and name of the file (with '') that stores the spectrum struct, or if already loaded, the variable name (without '')
% d        - [m] depth when the spectrum was generated, in orginal scale
% scaleUp  - [-] the length scale to apply
% sumA     - [-] desired linear amplitude sum, in target scale, set to 0 to leave it unmodified

% Li Ma (2013), lm808@ic.ac.uk, Imperial College London


if isstruct(filename)
    spec = filename;
else
    if strcmp(filename(end-2:end),'txt')
        fid = fopen(filename,'r');

        C = textscan(fid,'%f %f %f %f');

        spec.f = C{1};
        spec.a = C{2};
        spec.th = C{3};
        spec.ph = C{4};
        fclose(fid);

    warning('fSpecImport: Check txt file column orders!!');

    elseif strcmp(filename(end-2:end),'mat')
        load(filename)
    end
end

spec.a = spec.a * scaleUp;
spec.f = spec.f / sqrt(scaleUp);
spec.d = d * scaleUp;

spec.omega = 2*pi*spec.f;
spec.k = fDispersion(spec.d,1./spec.f);

if sumA~=0
    spec.a = spec.a/sum(spec.a)*sumA;
else
    disp('fSpecImport: amplitude unmodified.');
end
