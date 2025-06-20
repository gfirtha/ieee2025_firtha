function [position, Nout] = get_default_layout(Shape, Num_channels, R,z0, type)

switch Shape
    case 'Spherical'
        switch type
            case 'T'
                ssd = spherical_design('Design',Num_channels,'T',R);
                x = ssd.x0(:,1);
                y = ssd.x0(:,2);
                z = ssd.x0(:,3);
                Nout = size(ssd.x0,1);
        end
    case 'Circular'
        Nout = Num_channels;
    case 'Standard'
        Nout = Num_channels;
        switch Num_channels
            case 1
                [x,y,z] = sph2cart(0*pi/180 + pi/2,0,R);
            case 2
                [x,y,z] = sph2cart([30; -30]*pi/180 + pi/2,zeros(Num_channels,1),repmat(R,Num_channels,1));
            case 3
                [x,y,z] = sph2cart([30; -30; 0]*pi/180 + pi/2,zeros(Num_channels,1),repmat(R,Num_channels,1));
            case 4
                switch type
                    case 'Bformat'
                        [x,y,z] = sph2cart([45; -45; -135; -225]*pi/180 + pi/2,zeros(Num_channels,1),repmat(R,Num_channels,1));
                    case 'Dolby_stereo'
                        [x,y,z] = sph2cart([-30; 0; 30; 180]*pi/180 + pi/2,zeros(Num_channels,1),repmat(R,Num_channels,1));
                end
            case 5
                [x,y,z] = sph2cart([-30; 0; 30; 110; -110]*pi/180 + pi/2,zeros(Num_channels,1),repmat(R,Num_channels,1));
            case 6
                switch type
                    case 'Bformat'
                        [x1,y1,z1] = sph2cart([45; -45; -135; -225]*pi/180 + pi/4,zeros(4,1),repmat(R,4,1));
                        [x2,y2,z2] = sph2cart(zeros(2,1),[90;-90]*pi/180,repmat(R,2,1));
                        x = [x1;x2];
                        y = [y1;y2];
                        z = [z1;z2];
                end
            case 7
                switch type
                    case '7.0'
                        [x,y,z] = sph2cart([-30; 0; 30; 100; 140; -140; -100]*pi/180 + pi/2,zeros(Num_channels,1),repmat(R,Num_channels,1));
                end
            case 9
                switch type
                    case '5.0.4_Atmos'
                        [x1,y1,z1] = sph2cart([-30; 0; 30; 110; -110]*pi/180 + pi/2,zeros(5,1),repmat(R,5,1));
                        z_height = 1.2;
                        [x2,y2,z2] = sph2cart([-20; 20; 160; -160]*pi/180 + pi/2,zeros(4,1),repmat(R*0.5,4,1));
                        x = [x1;x2];
                        y = [y1;y2];
                        z = [z1;z2+z_height];
                    case '7.0.2_Atmos'
                        [x1,y1,z1] = sph2cart([-30; 0; 30; 100; 140; -140; -100]*pi/180 + pi/2,zeros(7,1),repmat(R,7,1));
                        z_height = 1.2;
                        [x2,y2,z2] = sph2cart([-45; 45]*pi/180 + pi/2,pi/180*45*ones(2,1),repmat(R*0.5,2,1));
                        x = [x1;x2];
                        y = [y1;y2];
                        z = [z1;z2+z_height];
                end

        end

end
position = [x,y,z+z0];

end

