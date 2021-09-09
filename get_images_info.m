function img_dir = get_images_info( s, os_cr, L, b, Nm, max_order )
%This function finds the information of the images resulting from a source
%in a rectangular room. 
%   Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London.
%       Input(s)
%-------------------------------------
%   s:  source location [x y x]
%   os_cr:   cart axes of oriented source [3x3] axis per row
%   L:  room dimension [x y z]
%   b:  walls reflection coefficient matrix [betha_x1 betha_y1 betha_z1;
%   betha_x2 betha_y2 betha_z2]
%   Nm:  maximum order of the wall images (N) [Nm_x Nm_y Nm_z]
%   Max_order: maximum order of reflection (-1 means the maximum possible order
%   based on Nm values)

%       Output(s)
%-------------------------------------
%   img_dir: image directory (cell) {order, p, m, Rpm, reflected_orientation, overall reflection coefficient }  -> Rpm (image location) and orientation are w/ respect to the origin 

img_dir=cell(0,6);  % cell {order, p, m, Rpm, reflected_orientation, overall reflection coefficient }

for px=0:1
    for py=0:1
        for pz=0:1
            for mx=-Nm(1):Nm(1)
                for my=-Nm(2):Nm(2)
                    for mz=-Nm(3):Nm(3)
                        m=[mx my mz];
                        p=[px py pz];
                        order=sum(abs(m-p))+sum(abs(m));
                        if ( (order <= max_order) || (max_order==-1) )
                            img_dir{end+1,1}=order;
                            img_dir{end,2}=p;
                            img_dir{end,3}=m;
                            img_dir{end,4}=(1-2*p).*s+2*m.*L;
                            %img_dir{end,5}=reflect_angle( s_or, p );
                            img_dir{end,5}=os_cr.*repmat((-1).^p,3,1);
                            m_minus_p=m-p;
                            img_dir{end,6}=prod(b(1,:).^abs(m_minus_p))*prod(b(2,:).^abs(m));
                        end
                        
                    end
                end
            end
        end
    end
end

end

