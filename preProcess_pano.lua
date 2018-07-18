-- sample code for saving gt images/edges/junctions/box parameters to .t7 file

require 'torch'
local matio = require 'matio'
matio.use_lua_strings = true -- read string
require 'image'
require 'paths'

im_path = './data/pano/img/' -- customize your gt data path
jc_path = './data/pano/junc/'
ed_path = './data/pano/edge/'
ln_path = './data/pano/line/'-- manhattan line
id_path = './data/pano/box_param/'

train_dir = paths.dir(im_path)

info_img_stack = torch.zeros(460, 3, 512, 1024) -- customize your gt image size
info_junc_stack = torch.zeros(460, 1, 512, 1024)
info_edge_stack = torch.zeros(460, 3, 512, 1024)
info_line_stack = torch.zeros(460, 3, 512, 1024)
info_id_stack = torch.zeros(460, 1024, 12)

cnt = 1
--fd = io.open('pano_rec_junc.txt', 'w')
for i = 1, #train_dir do

	dir = train_dir[i]
	if string.sub(dir, #dir - 2) == 'png' then
		print(dir)
		--fd:write(dir..'\n')

		-- img
		im = image.load(im_path..dir)
		info_img_stack[{{cnt},{},{},{}}] = im
		-- junc
		jc = image.load(jc_path..dir)
		info_junc_stack[{{cnt},{},{},{}}] = jc
                -- edg
                ed = image.load(ed_path..dir)
                info_edge_stack[{{cnt},{},{},{}}] = ed
                -- manhattan line
                ln = image.load(ln_path..dir)
                info_line_stack[{{cnt},{},{},{}}] = ln
                -- box param
                id = matio.load(id_path..string.sub(dir, 1, #dir - 4)..'.mat')
                id = id.box_n_all
                info_id_stack[{{cnt},{},{}}] = id

                cnt = cnt + 1
	end
end
print(cnt)

torch.save('./data/panoContext_img_train.t7', info_img_stack)
torch.save('./data/panoContext_edge_train.t7', info_edge_stack)
torch.save('./data/panoContext_line_train.t7', info_line_stack)
torch.save('./data/panoContext_cor_train.t7', info_junc_stack)
torch.save('./data/panoContext_box_train.t7', info_id_stack)

--fd:close()
