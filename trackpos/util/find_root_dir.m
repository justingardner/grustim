function rootdir = find_root_dir

    [ret, name] = system('hostname');
    if contains(name, 'Joshuas-MacBook-Air.local') ||contains(name, 'DN2t9suv.SUNet')  || contains(name, 'DNa80fac0.SUNet') 
        rootdir = '/Users/jryu/';
    elseif strncmp(name,'cointreau.stanford.edu',numel('cointreau.stanford.edu'))
        rootdir = '/Users/gru/';
    elseif ispc
        rootdir = winqueryreg('HKEY_CURRENT_USER',...
            ['Software\Microsoft\Windows\CurrentVersion\' ...
            'Explorer\Shell Folders'],'Personal');
    else
        rootdir = char(java.lang.System.getProperty('user.home'));
    end
    
end