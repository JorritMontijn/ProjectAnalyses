function [ChannelX,ChannelY] = MOL_GetChannelPosition(ProbeType,TipDepth)
%Script to get the depth of individual electrode sites given the probe
%configuration and the total depth of the probe. MOL2017
%Useful for CSD correction?

if isempty(TipDepth) || isnan(TipDepth)
    warning('No tipdepth was available from the excel overview: using 1000 as standard')
    TipDepth = 1000;
end

switch ProbeType
    
    case 'A1x32-Poly2-5mm-50s-177-CM32'
        ChannelX = NaN(32,1);
        ChannelY = NaN(32,1);
        
        TipToRow1 = 75;
        TipToRow2 = 100;
        
        InterSiteDistance = 50;
        
        Row1 = [16:-1:11 1:10];
        Row2 = [17:22 32:-1:23];
        
        Row1Depth = TipDepth - TipToRow1 - [0:15]*InterSiteDistance;
        Row2Depth = TipDepth - TipToRow2 - [0:15]*InterSiteDistance;

        ChannelY(Row1) = Row1Depth;
        ChannelY(Row2) = Row2Depth;

        ChannelX(Row1) = -43.3/2;
        ChannelX(Row2) = 43.3/2;

    case 'A4x8-5mm-100-200-177-CM32'
        ChannelX = NaN(32,1);
        ChannelY = NaN(32,1);
        
        TipToRow = 50;
        InterSiteDistance = 100;
        InterShankDistance = 200;

        Row1 = [1 8 2 7 3 6 4 5];
        Row2 = Row1 + 1*8;
        Row3 = Row1 + 2*8;
        Row4 = Row1 + 3*8;

        RowDepth = TipDepth - TipToRow - [0:7]*InterSiteDistance;

        ChannelY(Row1) = RowDepth;
        ChannelY(Row2) = RowDepth;        
        ChannelY(Row3) = RowDepth;
        ChannelY(Row4) = RowDepth;
        
        ChannelX(Row1) = 0*InterShankDistance;
        ChannelX(Row2) = 1*InterShankDistance;    
        ChannelX(Row3) = 2*InterShankDistance;
        ChannelX(Row4) = 3*InterShankDistance;
        
    case 'A1x64-Poly2-6mm-23s-160'
        
        ChannelX = NaN(64,1);
        ChannelY = NaN(64,1);
        
        TipToRow1 = 45;
        TipToRow2 = 22;
        
        InterSiteDistance = 46;
        
        Row1 = [32:-1:28 1:27];
        Row2 = [33:36 64:-1:37];
        
        Row1Depth = TipDepth - TipToRow1 - [0:31]*InterSiteDistance;
        Row2Depth = TipDepth - TipToRow2 - [0:31]*InterSiteDistance;

        ChannelY(Row1) = Row1Depth;
        ChannelY(Row2) = Row2Depth;

        ChannelX(Row1) = -30/2;
        ChannelX(Row2) = 30/2;
    otherwise
        error('Unknown probe type configuretion')
        
end
