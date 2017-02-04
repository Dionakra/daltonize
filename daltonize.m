%% Author: David de los Santos Boix <davsanboi@alum.us.es>
%% Keywords: daltonismo protanopia deuteranopia tritanopia
%% Maintainer: David de los Santos Boix <davsanboi@alum.us.es>
%% Created: 01/02/2017
%% Version: 1.0

%% usage: daltonize(IMAGE, [XYZ_TRANSFORMATION, LMS_TRANSFORMATION])
%%
%% Convert an ordinary image to protanopia, deuteranopia and tritanopia colorblind version.
%% 
%% XYZ_TRANSFORMATION list:
%%    adobe_rgb_d50     cie_rgb_d50       pal_secam_rgb_d50
%%    adobe_rgb_d65     cie_rgb_e         pal_secam_rgb_d65
%%    apple_rgb_d50     color_match_rgb   pro_photo_rgb
%%    apple_rgb_d65     don_rgb_4         smptec_rgb_d50
%%    best_rgb          eci_rgb           smptec_rgb_d65
%%    beta_rgb          ekta_space_ps5    srgb_d50
%%    bruce_rgb_d50     ntsc_rgb_c        srgb_d65
%%    bruce_rgb_d65     ntsc_rgb_d50      wide_gamut_rgb
%% 
%% LMS_TRANSFORMATION list:
%%    equal_energy
%%    equal_energy_d65
%%    ciecam97
%%    ciecam02
%%
%% Example:
%%  daltonize("beetlejuice.jpeg");
%%  daltonize("beetlejuice.jpeg", "wide_gamut_rgb");
%%  daltonize("beetlejuice.jpeg", "wide_gamut_rgb", "ciecam02");
function daltonize(image_path, xyz_transformation=0, lms_transformation=0)
  pkg load image;     % Image package needed
  
  % Checking if the image exists...
  if exist(image_path, 'file') == 0
    printf('The given image file does not exist\r\n');
    return;
  endif
  
  % Checking if the image is colorized...
  image = imread(image_path);
  sizeRGB = size(image);
  
  % If the image is in Black and White, then we cannot perform any colorblind operation on it
  if size(sizeRGB) != 3 || sizeRGB(3) != 3
    printf('The image should be colorized RGB\r\n');
    return;
  endif
  
  % Getting its values for Black/White and colorblind transformation
  BW = rgb2gray(image);
  RGB = double(image);
  
  % Choosing the transformation matrices...
  rgb2xyz = get_rgb2xyz_matrix(xyz_transformation);
  xyz2lms = get_xyz2lms_matrix(lms_transformation);

  % Now we calculate the necessary RGB to LMS and LMS to RGB matrices...
  rgb2lms = xyz2lms * rgb2xyz;
  lms2rgb = inv(rgb2lms);

  % Getting the LMS vectors configuration...
  lms_r = rgb2lms(1,:);
  lms_b = rgb2lms(3,:);
  lms_w = rgb2lms * ones(3, 1);

  % Applying the cross product to get the matrices...
  p0 = cross(lms_w, lms_b');
  p1 = cross(lms_w, lms_r');

  lms2lms_p = [  0.0000       0.0000,  0.0000;
               -p0(2)/p0(1),  1.0000,  0.0000;
               -p0(3)/p0(1),  0.0000,  1.0000]';
                    
  lms2lms_d = [ 1.0000, -p0(1)/p0(2),  0.0000;
                0.0000,       0.0000,  0.0000;
                0.0000, -p0(3)/p0(2),  1.0000]';

  lms2lms_t = [ 1.0000, 0.0000, -p1(1)/p1(3);
                0.0000, 1.0000, -p1(2)/p1(3);
                0.0000, 0.0000,      0.0000]';

  % Getting the colorblind images....
  tic
  RGB_p = get_colorblind_image(RGB, rgb2lms, lms2lms_p);
  RGB_d = get_colorblind_image(RGB, rgb2lms, lms2lms_d);
  RGB_t = get_colorblind_image(RGB, rgb2lms, lms2lms_t);
  toc
  
  % Formatting images to uint8...
  RGB_p = uint8(RGB_p);
  RGB_d = uint8(RGB_d);
  RGB_t = uint8(RGB_t);
 
  % Writing images to files...
  [dir, name, ext] = fileparts(image_path);
  imwrite(RGB_p,[name '_p' ext],'jpeg');
  imwrite(RGB_d,[name '_d' ext],'jpeg');
  imwrite(RGB_t,[name '_t' ext],'jpeg');
  imwrite(BW, [name '_bw' ext], 'jpeg');
endfunction

function RGB = get_colorblind_image(im, rgb2lms_matrix, colorblind_matrix)
  % First separate the RGB channels
  R = im(:,:,1);
  G = im(:,:,2);
  B = im(:,:,3);

  % Convert from RGB to LMS. The order is switched, so we have to transpose
  lms_image =  double([R(:), G(:), B(:)]) * rgb2lms_matrix';
  
  % Apply the colorblind matrix. Again, switched order, transposed matrix
  lms_cb = lms_image * colorblind_matrix';
  
  % Revert to RGB by inverting the rgb2lms_matrix
  RGB = lms_cb * inv(rgb2lms_matrix)' ;
  
  % Finally, reshape the image to its size
  RGB = reshape(RGB, size(im));  
endfunction

%http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
function rgb2xyz = get_rgb2xyz_matrix(type)
  switch type
    case 'adobe_rgb_d65'
      rgb2xyz = [ 0.5767309  0.1855540  0.1881852;
                  0.2973769  0.6273491  0.0752741;
                  0.0270343  0.0706872  0.9911085];
    case 'adobe_rgb_d50'
      rgb2xyz = [ 0.6097559  0.2052401  0.1492240;
                  0.3111242  0.6256560  0.0632197;
                  0.0194811  0.0608902  0.7448387];
    case 'apple_rgb_d65'
      rgb2xyz = [ 0.4497288  0.3162486  0.1844926;
                  0.2446525  0.6720283  0.0833192;
                  0.0251848  0.1411824  0.9224628];
    case 'apple_rgb_d50'
      rgb2xyz = [ 0.4755678  0.3396722  0.1489800;
                  0.2551812  0.6725693  0.0722496;
                  0.0184697  0.1133771  0.6933632];
    case 'best_rgb'
      rgb2xyz = [ 0.6326696  0.2045558  0.1269946;
                  0.2284569  0.7373523  0.0341908;
                  0.0000000  0.0095142  0.8156958];
    case 'beta_rgb'
      rgb2xyz = [ 0.6712537  0.1745834  0.1183829;
                  0.3032726  0.6637861  0.0329413;
                  0.0000000  0.0407010  0.7845090];
    case 'bruce_rgb_d65'
      rgb2xyz = [ 0.4674162  0.2944512  0.1886026;
                  0.2410115  0.6835475  0.0754410;
                  0.0219101  0.0736128  0.9933071];
    case 'bruce_rgb_d50'
      rgb2xyz = [ 0.4941816  0.3204834  0.1495550;
                  0.2521531  0.6844869  0.0633600;
                  0.0157886  0.0629304  0.7464909];
    case 'cie_rgb_e'
      rgb2xyz = [ 0.4887180  0.3106803  0.2006017;
                  0.1762044  0.8129847  0.0108109;
                  0.0000000  0.0102048  0.9897952];
    case 'cie_rgb_d50'
      rgb2xyz = [ 0.6343706  0.1852204  0.1446290;
                  0.3109496  0.5915984  0.0974520;
                 -0.0011817  0.0555518  0.7708399];
    case 'color_match_rgb'
      rgb2xyz = [ 0.5093439  0.3209071  0.1339691;
                  0.2748840  0.6581315  0.0669845;
                  0.0242545  0.1087821  0.6921735];
    case 'don_rgb_4'
      rgb2xyz = [ 0.6457711  0.1933511  0.1250978;
                  0.2783496  0.6879702  0.0336802;
                  0.0037113  0.0179861  0.8035125];
    case 'eci_rgb'
      rgb2xyz = [ 0.6502043  0.1780774  0.1359384;
                  0.3202499  0.6020711  0.0776791;
                  0.0000000  0.0678390  0.7573710];
    case 'ekta_space_ps5'
      rgb2xyz = [ 0.5938914  0.2729801  0.0973485;
                  0.2606286  0.7349465  0.0044249;
                  0.0000000  0.0419969  0.7832131];
    case 'ntsc_rgb_c'
      rgb2xyz = [ 0.6068909  0.1735011  0.2003480;
                  0.2989164  0.5865990  0.1144845;
                  0.0000000  0.0660957  1.1162243];
    case 'ntsc_rgb_d50'
      rgb2xyz = [ 0.6343706  0.1852204  0.1446290;
                  0.3109496  0.5915984  0.0974520;
                 -0.0011817  0.0555518  0.7708399];
    case 'pal_secam_rgb_d65'
      rgb2xyz = [ 0.4306190  0.3415419  0.1783091;
                  0.2220379  0.7066384  0.0713236;
                  0.0201853  0.1295504  0.9390944];
    case 'pal_secam_rgb_d50'
      rgb2xyz = [ 0.4552773  0.3675500  0.1413926;
                  0.2323025  0.7077956  0.0599019;
                  0.0145457  0.1049154  0.7057489];
    case 'pro_photo_rgb'
      rgb2xyz = [ 0.7976749  0.1351917  0.0313534;
                  0.2880402  0.7118741  0.0000857;
                  0.0000000  0.0000000  0.8252100];
    case 'smptec_rgb_d65'
      rgb2xyz = [ 0.3935891  0.3652497  0.1916313;
                  0.2124132  0.7010437  0.0865432;
                  0.0187423  0.1119313  0.9581563];
    case 'smptec_rgb_d50'
      rgb2xyz = [ 0.4163290  0.3931464  0.1547446;
                  0.2216999  0.7032549  0.0750452;
                  0.0136576  0.0913604  0.7201920];
    case 'srgb_d65'
      rgb2xyz = [ 0.4124564  0.3575761  0.1804375;
                  0.2126729  0.7151522  0.0721750;
                  0.0193339  0.1191920  0.9503041];
    case 'srgb_d50'
      rgb2xyz = [ 0.4360747  0.3850649  0.1430804;
                  0.2225045  0.7168786  0.0606169;
                  0.0139322  0.0971045  0.7141733];
    case 'wide_gamut_rgb'
      rgb2xyz = [ 0.7161046  0.1009296  0.1471858;
                  0.2581874  0.7249378  0.0168748;
                  0.0000000  0.0517813  0.7734287];
    otherwise
      % We choose wide_gamut_rgb as the default one
      type = 'wide_gamut_rgb';
      rgb2xyz = [ 0.7161046  0.1009296  0.1471858;
                  0.2581874  0.7249378  0.0168748;
                  0.0000000  0.0517813  0.7734287];
  endswitch
  
  printf('RGB to XYZ: %s\r\n', type);
endfunction

%https://en.wikipedia.org/wiki/LMS_color_space#XYZ_to_LMS
function xyz2lms = get_xyz2lms_matrix(type)
  switch type
    case 'equal_energy'
      xyz2lms = [ 0.38971  0.68898 -0.07868; 
                 -0.22981  1.18340  0.04641;
                  0.00000  0.00000  1.00000];
    case 'equal_energy_d65'
      xyz2lms = [ 0.4002  0.7076 -0.0808;
                 -0.2263  1.1653  0.0457;
                  0.0000  0.0000  0.9182];
    case 'ciecam97'
      xyz2lms = [ 0.8951  0.2664 -0.1614;
                 -0.7502  1.7135  0.0367;
                  0.0389 -0.0685  1.0296];
    case 'ciecam02'
      xyz2lms = [ 0.7328 0.4296 -0.1624;
                 -0.7036 1.6975  0.0061;
                  0.0030 0.0136  0.9834];
    otherwise
      % We choose CIECAM02 as the default one
      type = 'ciecam02';
      xyz2lms = [ 0.7328 0.4296 -0.1624;
                 -0.7036 1.6975  0.0061;
                  0.0030 0.0136  0.9834];
  endswitch
  printf('XYZ to LMS: %s\r\n', type);
endfunction