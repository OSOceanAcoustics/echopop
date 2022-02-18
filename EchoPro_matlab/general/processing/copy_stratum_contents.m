function  copy_stratum_contents(stratum0, stratum1, stratum2)
%% copy contents of stratum1 or averages of stratem1 and stratum2 to stratum0
%% created on 3/10/2020
%% 

global data

names = fieldnames(data.bio.strata);


switch nargin
    case 1
        fprintf('not enough input arguments, # of argument >= 2 ...\n')
        return
    case 2  % stratum1 to stratum0
        for i = 1:length(names)
            cmd = ['data.bio.strata(' num2str(stratum0) ').' char(names(i)) ' = data.bio.strata(' num2str(stratum1) ').' char(names(i)) ';'];
            eval(cmd)
        end
    case 3  % averages of stratum1 and stratum2 to stratum0
        for i = 1:length(names)
            cmd1 = ['data.bio.strata(' num2str(stratum0) ').' char(names(i)) ' = data.bio.strata(' num2str(stratum1) ').' char(names(i))];
            cmd2 = [' + data.bio.strata(' num2str(stratum2) ').' char(names(i)) ';'];
            cmd = [cmd1 cmd2];
            eval(cmd)
        end
    otherwise
        fprintf('not enough input arguments, # of argument <= 3 ...\n')
end

