-- driver

local cmd = torch.CmdLine()
cmd:text()
cmd:text('Script for training LayoutNet box parameter branch.')

cmd:option('-lr' , 1e-4, 'learning rate')
cmd:option('-epsilon' , 1e-6, 'adam, epsilon')
cmd:option('-batchSize' , 20, 'mini batch size')
cmd:option('-numPasses' , 1, 'number of passes')

cmd:text()
opt = cmd:parse(arg)

dofile('model_pano_box.lua')
dofile('train_pano_box.lua')
