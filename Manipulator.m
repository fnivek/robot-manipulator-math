classdef Manipulator
    %Manipulator Represents a manipulator such as a robotic arm
    %   Takes a DH table as its arguments, it can then produce all
    %   important attributes of the manipulator
    
    properties
        % Symbolic matrix that represents the DH table
        dh_table
        
        % Position and orientation of the end efector
        end_effector
        
        % Geometric jacobian
        geom_jacobian
        
        % Anylitical jacobian
        anl_jacobian
    end
    
    properties (Constant)
        THETA_INDEX = 1
        D_INDEX = 2
        A_INDEX = 3
        ALPHA_INDEX = 4
    end
    
    methods
        function obj = Manipulator(dh_table)
            %Constructor builds all the values based off dh_table
            %   dh_table is a symbolic matrix representing the dh_table
            %   should be ordered theta, d, a, alpha
            obj.dh_table = dh_table;
            obj.end_effector = GetTransform(obj, 0, size(dh_table, 1));
        end
        
        function T = GetTransform(obj, to, from)
            %Get the transform between to frame and from frame
            if to < 0 || to > size(obj.dh_table, 1) || from < 0 || from > size(obj.dh_table, 1)
                error('Bad frames')
            end
            
            T = eye(4);
            if to == from
                return
            end
            
            flip = 0;
            if to > from
                flip = 1;
                temp = to;
                to = from;
                from = temp;
            end
            
            for frame = [to + 1:from]
                dh_row = obj.dh_table(frame, :);
                ct = cos(dh_row(Manipulator.THETA_INDEX));
                ca = cos(dh_row(Manipulator.ALPHA_INDEX));
                st = sin(dh_row(Manipulator.THETA_INDEX));
                sa = sin(dh_row(Manipulator.ALPHA_INDEX));
                a = dh_row(Manipulator.A_INDEX);
                d = dh_row(Manipulator.D_INDEX);
                A = [ct, -st * ca, st * sa, a * ct; ...
                     st, ct * ca, -ct * sa, a * st; ...
                     0, sa, ca, d; ...
                     0, 0, 0, 1];
                T = T * A;
            end
            
            if flip == 1
                T = inv(T);
            end
        end
    end
    
end

