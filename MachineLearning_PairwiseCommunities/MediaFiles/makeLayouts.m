% This version of makeLayouts specifies the recipe for bifidobacterium media

addpath(genpath('/project/cometsfba/dimucci/comets-toolbox')) %% for use on scc1
addpath(genpath('/project/cometsfba/dimucci/cobratoolbox'))

cycles = 1000;
media = 0;
death = 0;
step = .1;
models = dir('*.mat');

for i = 1:length(models)
    x(i) = load(models(i).name);
end

for i = 1:length(models)
    y(i) = x(i).model;
    y(i).description = strrep(y(i).description,'.','_');
end

C = cell(1,length(y));
for i = 1:length(y)
    C{i} = y(i);
end
layout = createLayout(C{1:end});

% Modify parameters
layout.params.maxCycles =  cycles;
layout.params.deathRate = death;
layout.params.maxSpaceBiomass = 1;
layout.params.minSpaceBiomass = 0;
layout.params.timeStep = step;
layout.params.writeTotalBiomassLog = 'true';
layout.params.writeMediaLog = 'false';

% Initialize all media to 1e-6
layout.media_amt(1:end) = 0;

fid = fopen('compounds.m');
gid = fopen('ratios.m');
tline = fgetl(fid);
gline = fgetl(gid);
while ischar(tline)
        i = stridx(tline,layout.mets);
        if (isempty(i)) == 0
                layout.media_amt(i) = str2num(gline)*10;
                % maintain a static level of media for each element i.e. chemostat mode
                %layout.global_media_refresh(i,:) = str2num(gline);
        end
    tline = fgetl(fid);
    gline = fgetl(gid);
end

pop(1:length(y)) = 1e-9;
layout = setInitialPop(layout, '1x1', pop);



writeCometsLayout(layout,'.')
quit()
