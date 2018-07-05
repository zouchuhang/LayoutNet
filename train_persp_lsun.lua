-- require 'sys'
require 'image'

--local matio = require 'matio'
local image = require 'image'

sampleSize = opt.batchSize
numberOfPasses = opt.numPasses

function getBatch_val(data, sampsize, count)
    -- select batch
    inputMat = torch.zeros(sampsize, data.inp:size(2), data.inp:size(3), data.inp:size(4))
    gtMat = torch.zeros(sampsize, data.gt:size(2), data.gt:size(3), data.gt:size(4))
    for i = 1, sampsize do
        inputMat[{{i},{},{},{}}] = data.inp[{{count},{},{},{}}]
        data_gt = data.gt[{{count},{},{}, {}}]
        -- gaussian
        gtMat[{{i},{},{}, {}}] = data_gt
        count = count + 1
    end

    return inputMat, gtMat, count
end

function getBatch(data, sampsize, count, idx)
    -- select batch
    inputMat = torch.zeros(sampsize, data.inp:size(2), data.inp:size(3), data.inp:size(4))
    gtMat = torch.zeros(sampsize, data.gt:size(2), data.gt:size(3), data.gt:size(4))
    gtMsk = torch.zeros(sampsize, data.gt:size(2), data.gt:size(3), data.gt:size(4))

    for i = 1, sampsize do
        dir = data.tr_name[idx[count]]
        im = image.load(data.im_path..dir)
        data_inp = im
        ed = image.load(data.ed_path..dir)
        data_gt = ed
        -- data augmentation
        torch.seed() -- randomization
        -- flip
        local f_prob = torch.rand(1)
        if f_prob[1]>0.5 then
           data_inp = image.hflip(torch.reshape(data_inp, data.inp:size(2), data.inp:size(3), data.inp:size(4)))
           data_gt = image.hflip(torch.reshape(data_gt, data.gt:size(2), data.gt:size(3), data.gt:size(4)))  
           data_inp = torch.reshape(data_inp, 1, data.inp:size(2), data.inp:size(3), data.inp:size(4))
          data_gt = torch.reshape(data_gt, 1, data.gt:size(2), data.gt:size(3), data.gt:size(4))
	end
        -- gamma
        torch.seed() 
        local g_prob = torch.add(torch.mul(torch.rand(1),1.5), 0.5)
        data_inp = torch.pow(data_inp, g_prob[1])
        -- rotate
        torch.seed() 
        local r_prob = torch.rand(1)
        r_prob = r_prob[1]*1/9 - 1/18
        data_inp = image.rotate(torch.reshape(data_inp, data.inp:size(2), data.inp:size(3), data.inp:size(4)), r_prob, 'bilinear')
        data_gt = image.rotate(torch.reshape(data_gt, data.gt:size(2), data.gt:size(3), data.gt:size(4)), r_prob, 'bilinear')

        msk = data_gt:gt(0)
        inputMat[{{i},{},{},{}}] = data_inp
        gtMat[{{i},{},{}, {}}] = data_gt
        gtMsk[{{i},{},{}, {}}] = msk
        count = count + 1

        if count > tr_size then
            count = 1
            idx = torch.randperm(tr_size)
        end
    end

    return inputMat, gtMat, gtMsk, count, idx
end


function getValLoss()
    local valnumberOfPasses = torch.floor(pano_val.inp:size(1)/1)
    local loss = 0
    local valcount = 1
    --local out

    for i=1, valnumberOfPasses do

        --------------------- get mini-batch -----------------------
        inputMat, gtMat, valcount = getBatch_val(pano_val, 1, valcount)
        ------------------------------------------------------------      

        -- forward
        inputMat = inputMat:cuda()
        gtMat = gtMat:cuda()
        --print('forward')
        output = model.core:forward(inputMat)
       
        loss = model.criterion:forward(output, gtMat) + loss
        output = nil
        collectgarbage()
    end
    loss = loss / valnumberOfPasses

    return loss
end

-- do fwd/bwd and return loss, grad_params
function feval(x)
    if x ~= params then
        params:copy(x)
    end
    grad_params:zero()
    
    local loss = 0
    -- add for loop to increase mini-batch size
    for i=1, numberOfPasses do

        --------------------- get mini-batch -----------------------
        --inputMat, gtMat, gtMask = getBatch_rand(pano_tr, sampleSize)
        inputMat, gtMat, gtMsk, count, idx = getBatch(pano_tr, sampleSize, count, idx)
        ------------------------------------------------------------      

        -- forward
        inputMat = inputMat:cuda()
        gtMat = gtMat:cuda()
        
        output = model.core:forward(inputMat)

        --timer = torch.Timer()
        loss = model.criterion:forward(output, gtMat)-- + loss
        --print('Time elapsed : ' .. timer:time().real .. ' seconds')
        -- backward
        loss_d_1 = model.criterion:backward(output, gtMat)
        
        gtMsk = torch.mul(gtMsk, 4)
        gtMsk = gtMsk:cuda()
        gtMsk_w = torch.cmul(loss_d_1, gtMsk)
        loss_d = torch.add(gtMsk_w, loss_d_1)

        model.core:backward(inputMat, loss_d)

        output = nil
        loss_d = nil
        loss_d_1 = nil
        collectgarbage()
    end
    
    grad_params:div(numberOfPasses)
    
    -- clip gradient element-wise
    grad_params:clamp(-10, 10)
    
    
    return loss, grad_params
end

losses = {} 
vallosses = {}
local optim_state = {opt.lr, opt.epsilon}
--local optim_state = {learningRate = 1e-4, alpha = 0.95, epsilon = 1e-6}
local iterations = 8000
local minValLoss = 1/0
count = 1

idx = torch.randperm(pano_tr.inp:size(1))

for i = 1, iterations do
    
    
    model.core:training()
    local _, loss = optim.adam(feval, params, optim_state)
    --local _, loss = optim.rmsprop(feval, params, optim_state)

    print(string.format("update param, loss = %6.8f, gradnorm = %6.4e", loss[1], grad_params:clone():norm()))
    if i % 20 == 0 then
        print(string.format("iteration %4d, loss = %6.8f, gradnorm = %6.4e", i, loss[1], grad_params:norm()))
        model.core:evaluate()
        valLoss, output = getValLoss()
        vallosses[#vallosses + 1] = valLoss
        print(string.format("validation loss = %6.8f", valLoss))
        if minValLoss > valLoss then
            minValLoss = valLoss
            params_save = params:clone()
            nn.utils.recursiveType(params_save, 'torch.DoubleTensor')
            torch.save("./model/perspfull_edg_lsun.t7", params_save:double())
            print("------- Model Saved --------")
        end
        losses[#losses + 1] = loss[1]
    end
end
