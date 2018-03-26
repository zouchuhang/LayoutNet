#!/usr/bin/perl -w



####################

$input_dir = '/IUS/vmr7/dhoiem/context/images/keunho/ppm//';
$output_dir = '/IUS/vmr7/dhoiem/context/images/keunho/segments/';

$convertfrom = '.ppm';

@imagelist = glob $input_dir.'*'.$convertfrom;
print $input_dir.'*'.$convertfrom . "\n";

$count = 0;
foreach $image (@imagelist) {    
    $count = $count + 1;
    if ($count > 0) {
	$where = index($image, $convertfrom); 
	$where2 = index($image, '//');    
	$basename = substr($image, $where2+2, $where-$where2-2);
   
	print $basename .  "\n";
     	system('segment 0.8 100 100 ' . $image . ' ' . $output_dir . $basename . '.ppm');
    }

}
