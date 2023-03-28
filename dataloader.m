function[datasets]=dataloader(path)
files = dir(path);
dataset_names = [];
nb_files = 0;
for i=1:length(files)
   if length(files(i).name)>4 && isequal(files(i).name(end-3:end),'.mat')
        dataset_names = [dataset_names, files(i).name, ' ']; 
        nb_files = nb_files +1;
    end
end 
dataset_names = split(dataset_names);
datasets = struct;

for i=1:nb_files
    name = split(dataset_names(i),'.');
    first_split = split(name(1),'_');
    K_split = split(first_split(1),'=');
    K = str2num(cell2mat(K_split(2)));
    M_split = split(first_split(2),'=');
    M = str2num(cell2mat(M_split(2)));
    
    loaded_set = load(strcat('.\Datasets\',cell2mat(dataset_names(i))));
    datasets(i).data = loaded_set.Dataset;
    datasets(i).K = K;
    datasets(i).M = M;
    datasets(i).name = name(1);
end


end