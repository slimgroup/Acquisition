function jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, rseed, boatspeed, tfireint_min, tdelay, delayboat, figparams = nothing)

    #---------------------------------------------------------------------------------------------------------------------------------
    # jitter_airgunarraysJune outputs a structure array including the jittered acquisition parameters
    #
    # Use:
    #   jitacq = jitter_airgunarraysJune(ns, ds, dt, rndfactor, p, nboats, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig)
    
    # Input:
    #               ns - number of sources
    #               ds - source sampling interval
    #               dt - time sampling interval
    #        rndfactor - [rndfactor(1) rndfactor(2)]
    #                  - rndfactor(1): factor used to round the firing times (e.g., if firing times need to be rounded to the 
    #                    third decimal place then rndfactor(1) = 1000)
    #                  - rndfactor(2): factor used to round the jittered shot positions (e.g., if positions need to be rounded to the 
    #                    second decimal place then rndfactor(2) = 100)
    #                p - subsampling factor [NOTE: for data conventionally acquired at a source sampling of 50.0m 
    #                    (i.e., airguns fire every 20.0s), p = 2, 4, etc., when input data has ds = 25.0m, 12.5m, etc.]
    #           nboats - number of boats (or source vessels)
    #            rseed - random seed [NOTE: set two seeds for one boat (with two airgun arrays), and four seeds for two 
    #                    boats (with two airgun arrays each). Also, random seed should be chosen in a way such that the 
    #                    jittered shot numbers lie in the interval [1,ns].]
    #        boatspeed - speed of the boat (in meters/second)
    #     tfireint_min - minimal interval between jittered firing times (in seconds), which cannot be violated for a
    #                    pragmatic acquisition scenario (default is 10.0s)
    #           tdelay - time delay between airgun arrays on a boat (default is 10.0s)
    #        delayboat - time delay before the second boat starts firing (in seconds) [NOTE: when nboats = 1, delayboat = 0]
    #              fig - plot the acquisition scenario ('yes' or 'no')
    #        figparams - parameters to set figure properties, such as, font size, font weight, etc.
    
    
    # Author: Haneet Wason
    #         Seismic Laboratory for Imaging and Modeling
    #         Department of Earth, Ocean, and Atmospheric Sciences
    #         The University of British Columbia
    #         
    # Date: June, 2013
    
    # You may use this code only under the conditions and terms of the
    # license contained in the file LICENSE provided with this source
    # code. If you do not agree to these terms you may not use this
    # software.
    #---------------------------------------------------------------------------------------------------------------------------------
    
    # June 11, 2013
    # firing times rounded to third decimal place
    # sjit ---> NOT rounded (vs. rounded in jitter_airgunarrays) 
    
    # June 26, 2013
    # sjit ---> rounded to second decimal place
    
    # April 28, 2014
    # added paramter ---> figparams 
    
    
    # Jittered firing times and source positions for airgun array 1 on boat 1
    #rng(rseed(1));
    rng = MersenneTwister(2) #using Random
    dtfirejitb1arr1 = tfireint_min .+ rand!(rng,zeros(1,Int(round(ns/p))))*(2*tfireint_min);
    tfirejitb1arr1 = cumsum(dtfirejitb1arr1, dims = 2); 
    tfirejitb1arr1 = tfirejitb1arr1 .- tdelay;
    tfirejitb1arr1 = round.(rndfactor[1] .*tfirejitb1arr1) ./rndfactor[1];
    sjitb1arr1 = round.(rndfactor[2] .*boatspeed*tfirejitb1arr1) ./rndfactor[2];  
    gridsjitb1arr1IND = round.(sjitb1arr1./ds);
    
    println("Boat 1 - airgun array 1")
    println(string("Minimum interval between jittered firing times: " ,string(minimum(diff(tfirejitb1arr1,dims=2))) ," s"))
    println(string("Maximum interval between jittered firing times: ", string(maximum(diff(tfirejitb1arr1,dims=2))), " s"))
    if minimum(diff(tfirejitb1arr1,dims = 2)) < tfireint_min
       error("The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.")
    end
    
    println(string("First jittered shot number: ", string(minimum(gridsjitb1arr1IND))))
    println(string("Last jittered shot number: ", string(maximum(gridsjitb1arr1IND))))
    if minimum(gridsjitb1arr1IND) < 1 || ns < maximum(gridsjitb1arr1IND)
       error("The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.")
    end
    #fprintf('\n')
    
    
    # Jittered firing times and source positions for airgun array 2 on boat 1
    #rng(rseed(2));
    rng = MersenneTwister(2)
    dtfirejitb1arr2 = tfireint_min .+ rand!(rng, zeros(1,Int(round(ns/p)))) .*(2*tfireint_min);
    tfirejitb1arr2 = cumsum(dtfirejitb1arr2, dims = 2); 
    tfirejitb1arr2 = round.(rndfactor[1] .*tfirejitb1arr2) ./rndfactor[1];
    sjitb1arr2 = round.(rndfactor[2] .*boatspeed*tfirejitb1arr2) ./rndfactor[2];  
    gridsjitb1arr2IND = round.(sjitb1arr2 ./ds);
    
    println("Boat 1 - airgun array 2")
    println(string("Minimum interval between jittered firing times:",string(minimum(diff(tfirejitb1arr2,dims = 2))),"s"))
    println(string("Maximum interval between jittered firing times:", string(maximum(diff(tfirejitb1arr2,dims = 2))), "s"))
    if minimum(diff(tfirejitb1arr2,dims = 2)) < tfireint_min
       error("The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.")
    end
    
    println(string("First jittered shot number:", string(minimum(gridsjitb1arr2IND))))
    println(string("Last jittered shot number:", string(maximum(gridsjitb1arr2IND))))
    if minimum(gridsjitb1arr2IND) < 1 || ns < maximum(gridsjitb1arr2IND)
       error("The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.")
    end
    #fprintf('\n')
    
    
    # Jittered acquisition parameters for the second boat
    if nboats == 2
       # Jittered firing times and source positions for airgun array 1 on boat 2 
       rng = MersenneTwister(2)
       dtfirejitb2arr1 = tfireint_min .+ rand!(rng,zeros(1,Int(round(ns/p)))) .*(2*tfireint_min);
       tfirejitb2arr1 = cumsum(dtfirejitb2arr1,dims=2); 
       tfirejitb2arr1 = tfirejitb2arr1 .- tdelay;
       tfirejitb2arr1 = round.(rndfactor[1].*tfirejitb2arr1)/rndfactor[1];
       sjitb2arr1 = round.(rndfactor[2].*boatspeed.*tfirejitb2arr1) ./rndfactor[2];  
       gridsjitb2arr1IND = round.(sjitb2arr1./ds);
       tfirejitb2arr1 = delayboat .+ tfirejitb2arr1; 
       
       println("Boat 2 - airgun array 1")
       println(string("Minimum interval between jittered firing times:", string(minimum(diff(tfirejitb2arr1,dims=2))), "s"))
       println(string("Maximum interval between jittered firing times:", string(maximum(diff(tfirejitb2arr1,dims=2))),  "s"))
       if minimum(diff(tfirejitb2arr1,dims = 2)) < tfireint_min
          error("The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.")
       end
    
       println(string("First jittered shot number: ", string(minimum(gridsjitb2arr1IND))))
       println(string("Last jittered shot number: ", string(maximum(gridsjitb2arr1IND))))
       if minimum(gridsjitb2arr1IND) < 1 || ns < maximum(gridsjitb2arr1IND)
          error("The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.")
       end
       #fprintf('\n')
    
    
       # Jittered firing times and source positions for airgun array 2 on boat 2
       
       dtfirejitb2arr2 = tfireint_min .+ rand!(rng,zeros(1,Int(round(ns/p)))) .*(2 .*tfireint_min);
       tfirejitb2arr2 = cumsum(dtfirejitb2arr2,dims=2); 
       tfirejitb2arr2 = round.(rndfactor[1] .*tfirejitb2arr2)./rndfactor[1];
       sjitb2arr2 = round.(rndfactor[2] .*boatspeed .*tfirejitb2arr2)/rndfactor[2];     
       gridsjitb2arr2IND = round.(sjitb2arr2 ./ds);
       tfirejitb2arr2 = delayboat .+ tfirejitb2arr2; 
       
       println("Boat 2 - airgun array 2")
       println(string("Minimum interval between jittered firing times: ", string(minimum(diff(tfirejitb2arr2,dims=2))), " s"))
       println(string("Maximum interval between jittered firing times:", string(maximum(diff(tfirejitb2arr2,dims=2))), " s"))
       if minimum(diff(tfirejitb2arr2,dims = 2)) < tfireint_min
          error("The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.")
       end
    
       println(string("First jittered shot number: ", string(minimum(gridsjitb2arr2IND))))
       println(string("Last jittered shot number:", string(maximum(gridsjitb2arr2IND))))
       if minimum(gridsjitb2arr2IND) < 1 || ns < maximum(gridsjitb2arr2IND)
          error("The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.")
       end
       #fprintf('\n')
    end
    
    
    # Jittered firing time grid indices 
    if nboats == 1
       tfirejitgrid = 0 : dt : maximum([tfirejitb1arr1[end],tfirejitb1arr2[end]]);
       tfirejitb1arr1gridIND = zeros(1,length(tfirejitb1arr1)); 
       tfirejitb1arr2gridIND = zeros(1,length(tfirejitb1arr2)); 
       for j = 1:length(tfirejitb1arr1)
           ind_arr1 = findall(x-> x<=tfirejitb1arr1[j],tfirejitgrid);
           tfirejitb1arr1gridIND[j] = ind_arr1[end];
           ind_arr2 = findall(x-> x<=tfirejitb1arr2[j],tfirejitgrid);
           tfirejitb1arr2gridIND[j] = ind_arr2[end];
       end
       tshiftb1arr1 = tfirejitgrid[Int.(tfirejitb1arr1gridIND)] .- tfirejitb1arr1;
       tshiftb1arr2 = tfirejitgrid[Int.(tfirejitb1arr2gridIND)] .- tfirejitb1arr2;
    elseif nboats == 2
       gridend = [tfirejitb1arr1[end],tfirejitb1arr2[end],tfirejitb2arr1[end],tfirejitb2arr2[end]];
       tfirejitgrid = 0 : dt : maximum(gridend);
       tfirejitb1arr1gridIND = zeros(1,length(tfirejitb1arr1)); 
       tfirejitb1arr2gridIND = zeros(1,length(tfirejitb1arr2)); 
       tfirejitb2arr1gridIND = zeros(1,length(tfirejitb2arr1)); 
       tfirejitb2arr2gridIND = zeros(1,length(tfirejitb2arr2)); 
       for j = 1:length(tfirejitb1arr1)
           ind_b1arr1 = findall(x-> x<=tfirejitb1arr1[j],tfirejitgrid);
           tfirejitb1arr1gridIND[j] = ind_b1arr1[end];
           ind_b1arr2 = findall(x-> x<=tfirejitb1arr2[j],tfirejitgrid);
           tfirejitb1arr2gridIND[j] = ind_b1arr2[end];    
           ind_b2arr1 = findall(x-> x<=tfirejitb2arr1[j],tfirejitgrid);
           tfirejitb2arr1gridIND[j] = ind_b2arr1[end];
           ind_b2arr2 = findall(x-> x<=tfirejitb2arr2[j],tfirejitgrid);
           tfirejitb2arr2gridIND[j] = ind_b2arr2[end];
       end
       tshiftb1arr1 = tfirejitgrid[Int.(tfirejitb1arr1gridIND)] - tfirejitb1arr1;
       tshiftb1arr2 = tfirejitgrid[Int.(tfirejitb1arr2gridIND)] - tfirejitb1arr2;  
       tshiftb2arr1 = tfirejitgrid[Int.(tfirejitb2arr1gridIND)] - tfirejitb2arr1;  
       tshiftb2arr2 = tfirejitgrid[Int.(tfirejitb2arr2gridIND)] - tfirejitb2arr2;    
    end
    
    
    # Plot the acquisition scheme
    #if (@isdefined fig) || (@isdefined var) || fig == nothing
    #    fig = "no";
    #end
    
    #if cmp(fig, "yes") == 0
    #    figparams = Dict("fsize_label" => 16);
        #figparams.fsize_label = 16;
    #    figparams["fontname"] = "helvetica";
    #    figparams["fontweight"] = "demi";
    #    figparams["fsize_axis"] = 16;
       #if nboats == 1
       #   figure
       #   plot(sjitb1arr1, tfirejitb1arr1, 'o', sjitb1arr2, tfirejitb1arr2, 'b*'); 
       #   legend('Array 1', 'Array 2'); axis ij, axis('tight'); 
       #   xlabel('Source position (m)','FontSize',figparams.fsize_label,'FontName',figparams.fontname,'FontWeight',figparams.fontweight); 
       #   ylabel('Recording time (s)','FontSize',figparams.fsize_label, 'FontName',figparams.fontname,'FontWeight',figparams.fontweight);
       #   set(gca,'Fontsize',figparams.fsize_axis,'FontName',figparams.fontname,'FontWeight',figparams.fontweight,'TickDir','out')
    
       #elseif nboats == 2
       #   figure
       #   plot(sjitb1arr1, tfirejitb1arr1, 'o', sjitb1arr2, tfirejitb1arr2, 'b*', sjitb2arr1, tfirejitb2arr1, 'bo', sjitb2arr2, tfirejitb2arr2, 'b*'); 
       #   legend('Array 1', 'Array 2'); axis ij, axis('tight'); 
       #   xlabel('Source position (m)','FontSize',figparams.fsize_label,'FontName',figparams.fontname,'FontWeight',figparams.fontweight); 
       #   ylabel('Recording time (s)','FontSize',figparams.fsize_label, 'FontName',figparams.fontname,'FontWeight',figparams.fontweight);
       #   set(gca,'Fontsize',figparams.fsize_axis,'FontName',figparams.fontname,'FontWeight',figparams.fontweight,'TickDir','out')
    
       #end
    #end
    
    jitacq = Dict("tfirejitb1arr1" => tfirejitb1arr1);
    # Decide output variables
    if nboats == 1
       #jitacq["tfirejitb1arr1"]       = tfirejitb1arr1;
       jitacq["sjitb1arr1"]             = sjitb1arr1;
       jitacq["gridsjitb1arr1IND"]      = gridsjitb1arr1IND;
       jitacq["tfirejitb1arr1gridIND"]  = tfirejitb1arr1gridIND;
       jitacq["tshiftb1arr1"]           = tshiftb1arr1;
       jitacq["tfirejitb1arr2"]         = tfirejitb1arr2;
       jitacq["sjitb1arr2"]             = sjitb1arr2;
       jitacq["gridsjitb1arr2IND"]      = gridsjitb1arr2IND;
       jitacq["tfirejitb1arr2gridIND"]  = tfirejitb1arr2gridIND;
       jitacq["tshiftb1arr2"]           = tshiftb1arr2;
       jitacq["tfirejitgrid"]           = collect(tfirejitgrid');
    elseif nboats == 2
       #jitacq["tfirejitb1arr1"]         = tfirejitb1arr1;
       jitacq["sjitb1arr1"]            = sjitb1arr1;
       jitacq["gridsjitb1arr1IND"]      = gridsjitb1arr1IND;
       jitacq["tfirejitb1arr1gridIND"]  = tfirejitb1arr1gridIND;
       jitacq["tshiftb1arr1"]           = tshiftb1arr1;   
       jitacq["tfirejitb1arr2"]         = tfirejitb1arr2;
       jitacq["sjitb1arr2"]             = sjitb1arr2;
       jitacq["gridsjitb1arr2IND"]      = gridsjitb1arr2IND;
       jitacq["tfirejitb1arr2gridIND"]  = tfirejitb1arr2gridIND;
       jitacq["tshiftb1arr2"]           = tshiftb1arr2;   
       jitacq["tfirejitb2arr1"]         = tfirejitb2arr1;
       jitacq["sjitb2arr1"]             = sjitb2arr1;
       jitacq["gridsjitb2arr1IND"]      = gridsjitb2arr1IND;
       jitacq["tfirejitb2arr1gridIND"]  = tfirejitb2arr1gridIND;
       jitacq["tshiftb2arr1"]           = tshiftb2arr1;   
       jitacq["tfirejitb2arr2"]         = tfirejitb2arr2;
       jitacq["sjitb2arr2"]             = sjitb2arr2;
       jitacq["gridsjitb2arr2IND"]      = gridsjitb2arr2IND;
       jitacq["tfirejitb2arr2gridIND"]  = tfirejitb2arr2gridIND;
       jitacq["tshiftb2arr2"]           = tshiftb2arr2;   
       jitacq["tfirejitgrid"]           = collect(tfirejitgrid');
    end
    return jitacq
 # function end
end
    
