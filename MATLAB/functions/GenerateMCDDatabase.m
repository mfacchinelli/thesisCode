function GenerateMCDDatabase(ts,hs,D,folder)
%% Main File

%...Globals
global fullDatabase

%...Set up progress bar
i = 0;
wb = waitbar(0,'Setting up...','Name','MCD Download');

%...Number of files in total
N = length(ts)*length(hs);

%...Loop over files
for t = ts
    for h = hs
        %...Update waitbar
        i = i+1;
        wb = waitbar(i/N,wb,sprintf(['T: %i deg - H: %i km',newline,'%.1f %%'],t,h,i/N*100));

        %...Final file name
        dataFile = fullfile(folder,[num2str(t),'_',num2str(h),'.txt']);

        %...Download data
        commandString = '';
        for d = 1:ceil(D/4)
            rawFile = downloadMDCFiles(t,h,d);
            commandString = horzcat(commandString,...
                sprintf('cat "%s" >> "%s";',rawFile,dataFile));
        end

        %...Append files
        system(['touch "',dataFile,'";',commandString]);
    end
end

%...Close progress bar
close(wb);
    
%% Supporting Functions

    function outFile = downloadMDCFiles(t,h,mode)

        switch mode
            case 1
                vars = 'var1=rho&var2=p&var3=t&var4=wind';
                saveFile = fullfile(folder,'first.txt');
            case 2
                if fullDatabase
                    vars = 'var1=w&var2=cp&var3=co2&var4=ar';
                else
                    vars = 'var1=w&var2=cp&var3=none&var4=none';
                end
                saveFile = fullfile(folder,'second.txt');
            case 3
                vars = 'var1=n2&var2=co&var3=o3&var4=o';
                saveFile = fullfile(folder,'third.txt');
            case 4
                vars = 'var1=o2&var2=hydro&var3=hydro2&var4=he';
                saveFile = fullfile(folder,'forth.txt');
            case 5
                vars = 'var1=h2ovap&var2=none&var3=none&var4=none';
                saveFile = fullfile(folder,'fifth.txt');
        end

        %...File names
        tempFile = fullfile(folder,'temporary.txt');

        %...Get URL for data
        url = ['http://www-mars.lmd.jussieu.fr/mcd_python/cgi-bin/mcdcgi.py?datekeyhtml=1&ls=',num2str(t),...
               '&localtime=0.&year=2018&month=3&day=13&hours=21&minutes=58&seconds=30&julian=2458191.415625&',...
               'martianyear=34&martianmonth=5&sol=304&latitude=all&longitude=all&altitude=',num2str(h*1e3),...
               '.&zkey=2&isfixedlt=off&dust=1&hrkey=1&zonmean=off&',vars,'&dpi=eps&islog=off&',...
               'colorm=jet&minval=&maxval=&proj=cyl&plat=&plon=&trans=&iswind=off&latpoint=&lonpoint='];

        %...Download data
        fileName = websave(tempFile,url,weboptions('Timeout',60));

        %...Open and read downloaded file
        fileID = fopen(fileName,'r');
        data = textscan(fileID,'%s');
        fclose(fileID);

        %...Extract txt file name
        filePath = data{1}{8}(10:end-7);

        %...Download new file
        url = 'http://www-mars.lmd.jussieu.fr/mcd_python/';
        outFile = websave(saveFile,[url,filePath]);

    end

end