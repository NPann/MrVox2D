clear A P
for a = 1:numel(Seq.RF.exc.phainc)
%     close all
    [A(a,:) P(a,:)] = VoxelSim2D_do_one('ModelPar.txt','SeqPar_TrueFISP.txt',1,a);
end

[A P] = VoxelSim2D_do_one('ModelPar.txt','SeqPar_Rand.txt',1,1);
 [A2(a,:) P2(a,:)] = VoxelSim2D_do_one('ModelPar.txt','SeqPar_TrueFISP.txt',1,a);