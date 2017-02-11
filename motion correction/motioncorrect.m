function motioncorrect

clear
addpath('/media/data/DATA1/Code/NoRMCorre-master');

dirpath = uigetdir; % do "uigetdir('your/home/directory/') if you want to make this better on your personal account
fs = strfind(dirpath, '/');
foldername = dirpath(fs(end)+1:end); clear fs

donotproceed = 1;

%

while donotproceed == 1
    prompt = ['Are you sure you want to proceed? The TIF sequence within ' dirpath ' will be processed. [y]/n: '];
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

    p = dir(dirpath);
    
    tic
    filenames={p.name};
    tif_files={};
    for j=1:size(filenames,2)
        if endsWith(filenames{j},'ome.tif')
            tif_files=[tif_files; filenames{j}];
        end
    end
    
    if ~isempty(tif_files)
        
        tmp=loadtiff([dirpath '/' tif_files{1}]);
        tif_stack = uint16(zeros(size(tmp, 1), size(tmp, 2), size(tif_files,1)));
        for j=1:size(tif_files,1)
            tif = loadtiff([dirpath '/' tif_files{j}]);
            
            tif_stack(:,:,j) = tif;
            
            %         tif_stack=cat(3, tif_stack, tif);
            delete([dirpath '/' tif_files{j}])
        end
        
        saveastiff(tif_stack, [dirpath '.tif']);
        fprintf(['The tif series ' foldername ' have been converted to a tif stack. Elapsed time : %.3f s.\n'], toc);
        fprintf('Starting NoRMCorre\n');
        %
    else
        fprintf(['There is no tif sequence in ' foldername '... Moving on...\n']);
    end
    
    
    if ~(0 == exist([dirpath '.tif'])) && (0 == exist([dirpath '_normcorre.tif']))
        tic
        [tif_stack_corrected, shifts] = normcorre_kheirbek([dirpath '.tif']);
        saveastiff(tif_stack_corrected, [dirpath '_normcorre.tif']);
        save([dirpath '_normcorre_shifts'], 'shifts');
        
        fprintf([ foldername ' has been motion corrected and saved as ' foldername '_normcorre.tif. Elapsed time : %.3f s.\n'], toc);
    elseif ~(0 == exist([dirpath '.tif'])) && ~(0 == exist([dirpath '_normcorre.tif']))
        
        fprintf(['TIF stack ' foldername '.tif and ' foldername '_normcorre.tif both exist. No work to do here...\n']);
    else
        fprintf(['TIF stack ' foldername '.tif does not exist... Moving on...\n']);
    end
    
elseif donotproceed == 3
    fprintf('NoRMCorre terminated\n');
end

end
