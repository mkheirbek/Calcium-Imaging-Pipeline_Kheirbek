function motioncorrect_batch

clear
addpath('/media/data/DATA1/Code/NoRMCorre-master');

dirpath = uigetdir; % do "uigetdir('your/home/directory/') if you want to make this better on your personal account

donotproceed = 1;

while donotproceed == 1
    prompt = ['Are you sure you want to proceed? ALL EXPERIMENTAL FOLDERS within ' dirpath ' will be processed. [y]/n: '];
    areyousure = input(prompt, 's');
    if areyousure == 'y'
        donotproceed = 2;
    elseif areyousure == 'n'
        donotproceed = 3;
    else
        donotproceed = 1;
    end
end

%
if donotproceed == 2
    
    fprintf('NoRMCorre is continuing...\n');
    
    %
    p = dir(dirpath);
    folderpaths = {p(3:end).name};
    folderpaths = folderpaths([p(3:end).isdir]);
    
    
    %
    
    for i = 1:size(folderpaths, 2)
        % Sebi's tif series loading, saving, and deleting code. Save the tif
        % stack with the same name as the folder
        tic
        subdirpath=dir([dirpath '/' folderpaths{i}]);
        filenames={subdirpath.name};
        
        tif_files={};
        for j=1:size(filenames,2)
            if endsWith(filenames{j},'ome.tif')
                tif_files=[tif_files; filenames{j}];
            end
        end
        
        if ~isempty(tif_files)
            
            tmp=loadtiff([dirpath '/' folderpaths{i} '/' tif_files{1}]);
            tif_stack = uint16(zeros(size(tmp, 1), size(tmp, 2), size(tif_files,1)));
            for j=1:size(tif_files,1)
                tif = loadtiff([dirpath '/' folderpaths{i} '/' tif_files{j}]);
                
                tif_stack(:,:,j) = tif;
                
                %         tif_stack=cat(3, tif_stack, tif);
                delete([dirpath '/' folderpaths{i} '/' tif_files{j}])
            end
            
            saveastiff(tif_stack, [dirpath '/' folderpaths{i} '/' folderpaths{i} '.tif']);
            fprintf(['The tif series ' folderpaths{i} ' have been converted to a tif stack. Elapsed time : %.3f s.\n'], toc);
            fprintf('Starting NoRMCorre\n');
            %
        else
            fprintf(['There is no tif sequence in ' folderpaths{i} '... Moving on...\n']);
        end
        
        
        if ~(0 == exist([dirpath '/' folderpaths{i} '/' folderpaths{i} '.tif'])) && (0 == exist([dirpath '/' folderpaths{i} '/' folderpaths{i} '_normcorre.tif']))
            tic
            [tif_stack_corrected, shifts] = normcorre_kheirbek([dirpath '/' folderpaths{i} '/' folderpaths{i} '.tif']);
            saveastiff(tif_stack_corrected, [dirpath '/' folderpaths{i} '/' folderpaths{i} '_normcorre.tif']);
            save([dirpath '/' folderpaths{i} '/' folderpaths{i} '_normcorre_shifts'], 'shifts');
            
            fprintf([ folderpaths{i} ' has been motion corrected and saved as ' folderpaths{i} '_normcorre.tif. Elapsed time : %.3f s.\n'], toc);
        elseif ~(0 == exist([dirpath '/' folderpaths{i} '/' folderpaths{i} '.tif'])) && ~(0 == exist([dirpath '/' folderpaths{i} '/' folderpaths{i} '_normcorre.tif']))
            
            fprintf(['TIF stack ' folderpaths{i} '.tif and ' folderpaths{i} '_normcorre.tif both exist. No work to do here...\n']);
        else
            fprintf(['TIF stack ' folderpaths{i} '.tif does not exist... Moving on...\n']);
        end
        
        
        
    end
    
    fprintf('All tif series in folder have been motion corrected through NoRMCorre\n');
    
elseif donotproceed == 3
    fprintf('NoRMCorre terminated\n');
end

end
