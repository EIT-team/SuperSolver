function write_electrode_nodes(elec,filename)

    for i=1:size(elec,1)
        elecnodes = unique(elec(i,:));
        if elecnodes(1) == 0
            elecnodes = elecnodes(2:end);
        end
        if i==1
            dlmwrite(filename,elecnodes);
        else
            dlmwrite(filename,elecnodes,'-append');
        end
    end