%{
 	Evaluates the bacterial cell envelope parameters
 	Parameters: 
        - @width, initial width in micrometers.
        - @length, initial length in micrometers.
        - @max, max lenght in micrometers.

    Note: 
        Use initial width == length to simulate spherical cells.
%}
function S = BctSimulateCellGrowth(length, width, max)
    w= width;               % initial width in micrometers.
	l= length;              % initial length in micrometers.
    avgd= 1000.0;           % Average cell density in kg/m^3.
    stdd= 4.0;              % Standard deviation of cell density in kg/m^3.
    
    
    delta= (max - l)/100;
    
    S= zeros(1000,7);
    i= 1;
    
    l1= 0;
    while l < max
        v= BctVolume(l,w);
        s= BctArea(l,w);
        mass= normrnd(avgd,stdd) * BctVolume(l*micrometer,w*micrometer);
        
        S(i,:) = [ i, l, v, s, (s/v), mass/femtogram, 0]; 
        
        l0= l;
        w0= w;
        l= l + delta;
        if width == length
            w= l;
        end   
        
        w1= w;
        if( i > 1 )
            m1= mass;
            l1= (w1^3 * m1 - w0^3 * m1 + 3 * w0^2 * l0 * m1) / (3 * m0 * w1^2);
        end
        m0= mass;
        
        S(i,7) = l1;     
        
        i= i + 1;
    end
    S= S(1:i-1,:);
    BctPlotCellData(S); 
end

