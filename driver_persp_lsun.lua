-- driver

local cmd = torch.CmdLine()

cmd:option('-lr' , 1e-4, 'learning rate')
cmd:option('-epsilon' , 1e-6, 'adam, epsilon')
cmd:option('-batchSize' , 10, 'mini batch size')
cmd:option('-numPasses' , 1, 'number of passes')

cmd:text()
opt = cmd:parse(arg)

dofile('model_persp_lsun.lua')
dofile('train_persp_lsun.lua')
