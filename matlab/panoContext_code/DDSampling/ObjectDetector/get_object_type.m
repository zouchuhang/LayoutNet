function [ output ] = get_object_type( input )
%GET_OBJECT_TYPE Summary of this function goes here
%   Detailed explanation goes here
% try 
%     I = isa(input, 'cell') & isa(input{1}, 'char');
% catch
%     fprintf('');
% end
if isa(input, 'cell') && isa(input{1}, 'char')
    output = zeros(length(input),1);
    for i = 1:length(input)
        switch input{i}
            case {'painting', 'picture', 'blanket'}
                output(i) = 1;
            case {'bedside','beside'}
                output(i) = 2;
            case {'door','doorway'}
                output(i) = 3;
            case {'mirror'}
                output(i) = 4;
            case {'cabinet'}
                output(i) = 5;
            case {'tv'}
                output(i) = 6;
            case {'desk','table'}
                output(i) = 7;
            case {'window'}
                output(i) = 8;
            case {'bed'}
                output(i) = 9;
            case {'chair'}
                output(i) = 10;
            case {'wardrobe'}
                output(i) = 11;
            case {'sofa'}
                output(i) = 12;
            case {'cushion'}
                output(i) = 13;
            case {'tv stand'}
                output(i) = 14;
            case {'bench'}
                output(i) = 15;
            case {'windowsill'}
                output(i) = 16;
            case {'trash can','basket'}
                output(i) = 17;
            case {'refrigerator'}
                output(i) = 18;
            case {'column','pillar'}
                output(i) = 19;
            case {'aircon'}
                output(i) = 20;
            case {'stair'}
                output(i) = 21;
            case {'washtub'}
                output(i) = 22;
            case {'bedhead'}
                output(i) = 23;
            case {'microwave oven'}   
                output(i) = 24;
            case {'heater'}
                output(i) = 25;
            case {'box','cuboid'}
                output(i) = 26;
            case {'shelf'}
                output(i) = 27;
            case {'room'}
                output(i) = 29;
            otherwise
                output(i) = 28;
        end
    end
    
elseif isa(input, 'numeric')
    output = cell(length(input),1);
    for i = 1:length(input)
        switch input(i)
            case 0
                output{i} = 'background';
            case 1
                output{i} = 'painting';
            case 2
                output{i} = 'bedside';
            case 3
                output{i} = 'door';
            case 4
                output{i} = 'mirror';
            case 5
                output{i} = 'cabinet';
            case 6
                output{i} = 'tv';
            case 7
                output{i} = 'desk';
            case 8
                output{i} = 'window';
            case 9
                output{i} = 'bed';
            case 10
                output{i} = 'chair';
            case 11
                output{i} = 'wardrobe';
            case 12
                output{i} = 'sofa';
            case 13
                output{i} = 'cushion';
            case 14
                output{i} = 'tv stand';
            case 15
                output{i} = 'bench';
            case 16
                output{i} = 'windowsill';
            case 17
                output{i} = 'trash can';
            case 18
                output{i} = 'refrigerator';
            case 19
                output{i} = 'column';
            case 20
                output{i} = 'aircon';
            case 21
                output{i} = 'stair';
            case 22
                output{i} = 'washtub';
            case 23
                output{i} = 'bedhead';
            case 24 
                output{i} = 'microwave oven';
            case 25
                output{i} = 'heater';
            case 26
                output{i} = 'cuboid';
            case 27
                output{i} = 'shelf';
            case 29
                output{i} = 'room';
            otherwise
                output{i} = 'other';
        end
    end
    
end

end

