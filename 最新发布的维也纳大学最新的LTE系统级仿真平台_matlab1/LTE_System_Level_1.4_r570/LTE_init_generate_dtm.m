function digital_terrain_model = LTE_init_generate_dtm(dtm_folder_,dtm_file_name,dtm_hdr_file_name,enable_plotting)
    % Martin Taranetz, INTHFT 2009
    
    % Function returns elevation map and description data with coordinates,
    % resolution and matrix size
    
    %% Initial Parameters
    enable_plotting = false;
    
    %% Read out the DTM .bil file and generate elevation map
    dtmfile_    = fopen(fullfile(dtm_folder_, dtm_file_name), 'r');
    dtmhdrfile_ = fopen(fullfile(dtm_folder_, dtm_hdr_file_name), 'r');
    
    dtm_hdr = textscan(dtmhdrfile_, '%s %f');
    
    dtm.description.NWxmap = dtm_hdr{2}(1);
    dtm.description.NWymap = dtm_hdr{2}(2);
    dtm.description.xdim   = dtm_hdr{2}(3);
    dtm.description.ydim   = dtm_hdr{2}(4);
    dtm.description.ncols  = dtm_hdr{2}(5);
    dtm.description.nrows  = dtm_hdr{2}(6);
    dtm.description.nbits  = dtm_hdr{2}(7);
    dtm.description.nbands = dtm_hdr{2}(8);
    
    dtm.description.SWxmap = dtm.description.NWxmap;                                                      % Lower-leftmost corner (x) -> SW
    dtm.description.SWymap = dtm.description.NWymap - dtm.description.ydim*(dtm.description.nrows-1/100); % Lower-leftmost corner (y) -> SW
    
    dtm.description.NExmap = dtm.description.NWxmap + dtm.description.xdim*(dtm.description.ncols-1/100);
    dtm.description.NEymap = dtm.description.NWymap;
    
    dtm.description.SExmap = dtm.description.NExmap;
    dtm.description.SEymap = dtm.description.SWymap;
    
    dtm.description.roi_x = [dtm.description.SWxmap dtm.description.SExmap];
    dtm.description.roi_y = [dtm.description.SWymap dtm.description.NWymap];
    
    if dtm.description.nbits == 16
        % Read in Matlabs columnwise manner and transpose
        dtm_map_t = fread(dtmfile_,[dtm.description.ncols, dtm.description.nrows],'int16');
    end
    dtm.data  = flipud(dtm_map_t');
    
    if enable_plotting
        figure;
        imagesc(dtm.description.roi_x,dtm.description.roi_y,dtm.data);
        hold on;
        set(gca,'YDir','normal');
        xlabel('x pos [m]');
        ylabel('y pos [m]');
        colorbar;
        title('Elevation map (m)');
    end
    
    digital_terrain_model = dtm;
    
    fclose(dtmfile_);
    fclose(dtmhdrfile_);
end