-- driver

local cmd = torch.CmdLine()
cmd:text()
cmd:text('Script for training sequence model.')

cmd:option('-lr' , 1e-4, 'learning rate')
cmd:option('-epsilon' , 1e-6, 'adam, epsilon')
cmd:option('-batchSize' , 1, 'mini batch size')
cmd:option('-numPasses' , 10, 'number of passes')
cmd:text()
opt = cmd:parse(arg)

dofile('model_pano_joint_lsun.lua')
dofile('train_pano_joint_lsun.lua')
