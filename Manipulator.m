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
        % DH table indices
        THETA_INDEX = 1
        D_INDEX = 2
        A_INDEX = 3
        ALPHA_INDEX = 4
        JOINT_TYPE_INDEX = 5

        % Joint types
        PRISMATIC = 0
        REVOLUTE = 1
        STATIC = 2

        % Type of Analytical Jacobian
        ZYZ_BODY_FIXED = 0
        ZYZ_SPACE_FIXED = 1
    end
    
    methods
        function obj = Manipulator(dh_table)
            %Constructor builds all the values based off dh_table
            %   dh_table is a symbolic matrix representing the dh_table
            %   should be ordered theta, d, a, alpha
            obj.dh_table = dh_table;
            obj.end_effector = GetTransform(obj, 0, size(dh_table, 1));
            obj.geom_jacobian = GetGeomJacobian(obj);
            obj.anl_jacobian = GetAnalyticalJacobian(obj, Manipulator.ZYZ_BODY_FIXED);
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

        function J = GetGeomJacobian(obj)
            J = sym(zeros(6, size(obj.dh_table, 1)));

            % Iterate over all frames 
            for frame = [1:size(obj.dh_table, 1)]
                dh_row = obj.dh_table(frame, :);
                joint_type = dh_row(Manipulator.JOINT_TYPE_INDEX);
                % Transform to the i-1 frame
                T_frame = GetTransform(obj, 0, frame - 1);
                % Z axis of the i - 1 frame
                z = T_frame(1:3, 3);
                end_effector_origin = obj.end_effector(1:3, 4);
                if joint_type == Manipulator.PRISMATIC
                    % Add prismatic column
                    % Linear velocity: z-axis_(i-1)
                    J(1:3, frame) = z;
                    % Angular velocity: 0
                elseif joint_type == Manipulator.REVOLUTE
                    % Add revolute column
                    % Linear velocity: z-axis(i-1) X (o_(n) - o_(i-1))
                    frame_origin = T_frame(1:3, 4);
                    r = end_effector_origin - frame_origin;
                    J(1:3, frame) = cross(z, r);

                    % Angular Velocity: z-axis(i-1)
                    J(4:6, frame) = z;
                end 
                % Else add 0 column for static manipulator
            end
        end

        function Ja = GetAnalyticalJacobian(obj, form)
            % Produce the analytical jacobian of form form
            psi = sym('psi');
            theta = sym('theta');
            T = zeros(3);
            if form == Manipulator.ZYZ_SPACE_FIXED
                T = [cos(psi)*sin(theta), -sin(psi), 0; sin(psi)*sin(theta), cos(psi), 0; cos(theta) 0 1];
            elseif form == Manipulator.ZYZ_BODY_FIXED
                T = [0, -sin(psi), cos(psi)*sin(theta); 0, cos(psi), sin(psi)*sin(theta); 1 0 cos(theta)];
            end
            Ja = [eye(3), zeros(3); zeros(3), inv(T)] * obj.geom_jacobian;
        end
    end
    
end