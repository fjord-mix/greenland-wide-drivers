figure; hold on
for i1=1:size(ensemble,2)
    for i2=1:size(ensemble,3)
    for i3=1:size(ensemble,4)
    for i4=1:size(ensemble,5) % we skip 'm' because it is being used for the metadata structure elsewhere in the code
    for i5=1:size(ensemble,6)
        plot(tf_box_comp(:,i1,i2,i3,i4,i5),-res_obs(i_fjord).zf);
    end
    end
    end
    end
end