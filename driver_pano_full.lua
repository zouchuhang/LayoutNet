-- driver

local cmd = torch.CmdLine()
cmd:text()
cmd:text('Script for training LayoutNet.')

cmd:option('-lr' , 1e-4, 'learning rate')
cmd:option('-epsilon' , 1e-6, 'adam, epsilon')
cmd:option('-batchSize' , 1, 'mini batch size')
cmd:option('-numPasses' , 20, 'number of passes, works as batch size to avoid out of memory')

cmd:text()
opt = cmd:parse(arg)

dofile('model_pano_full.lua')
dofile('train_pano_full.lua')
