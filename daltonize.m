%% Author: David de los Santos Boix <davsanboi@alum.us.es>
%% Keywords: daltonismo protanopia deuteranopia tritanopia
%% Maintainer: David de los Santos Boix <davsanboi@alum.us.es>
%% Created: 01/02/2017
%% Version: 1.0

%% usage: daltonize(IMAGE, [XYZ_TRANSFORMATION, LMS_TRANSFORMATION])
%%
%% Convert an ordinary image to a given colorblind type.
%% The user can customize the transformation matrices in order to 
%% play with them and view the different results applying different
%% transformation matrices.
%%
%% Example:
%%  daltonize("beetlejuice.jpeg");
%%  daltonize("beetlejuice.jpeg", "wide_gamut_rgb");
%%  daltonize("beetlejuice.jpeg", "wide_gamut_rgb", "ciecam02");

function daltonize(image_path, xyz_transformation, lms_transformation)
  % First of all we load the 'image' package to prevent lack of libraries
  pkg load image;

  % Here we load the chosen transformations by the user
  %   If the user doesn't provide any or mistaken ones, we chose the better for them
  if nargin == 1        
    rgb2xyz = get_rgb2xyz_matrix("wide_gamut_rgb");
    xyz2lms = get_xyz2lms_matrix("ciecam02");
    
  elseif nargin == 2
    xyz2lms = get_xyz2lms_matrix("ciecam02");
    rgb2xyz = get_rgb2xyz_matrix(xyz_transformation);
    
    if rgb2xyz == 0
      rgb2xyz = get_rgb2xyz_matrix("wide_gamut_rgb");
    endif
    
  elseif nargin == 3
    rgb2xyz = get_rgb2xyz_matrix(xyz_transformation);
    if rgb2xyz == 0
      rgb2xyz = get_rgb2xyz_matrix("wide_gamut_rgb");
    endif
    
    xyz2lms = get_xyz2lms_matrix(lms_transformation);
    if xyz2lms == 0
      xyz2lms = get_xyz2lms_matrix("ciecam02");
    endif
  else
    print("Incorrect number of input parameters");
    exit(-1);
  endif
  
  
  % Now we check if the given file exists
  if exist(image_path, "file") == 0
    print("The given image file doesn't exist");
    exit(-1);
  endif
        
  % This is the original rgb2lms transform matrix, which calculation I personally don't know where it comes from
  %   So, we will use the mathematic method, described in whe Android SDK, just here:
  %     https://android.googlesource.com/platform/frameworks/native/+/14e8b01/services/surfaceflinger/Effects/Daltonizer.cpp
  %   rgb2lms = [17.8824 43.5161 4.11935; 3.45565 27.1554 3.86714; 0.0299566 0.184309 1.46709]

  % Now we calculate the RGB2LMS transform matrix. The process is the following:
  %   XYZ = M * RGB
  %   LMS = M' * XYZ
  % So, finally, the transform matrix is:
  %   LMS = M' * M * RGB
  % So, in order to hurry the things up, we calculate the M'* M matrix just now and its inverse
  rgb2lms = xyz2lms * rgb2xyz;
  lms2rgb = inv(rgb2lms);

  % Finally, we load the image and get its RGB and BW codification
  image = imread(image_path);
  BW = rgb2gray(image);
  RGB = double(image);
  sizeRGB = size(RGB);
  
  % After loading the image, we need to split the image in its channels
  % It is due to the lack of perception of the blindcolor people to certain channels
  % This is explained in the following link:
  %   https://en.wikipedia.org/wiki/LMS_color_space
  %   https://github.com/joergdietrich/daltonize/blob/master/doc/project_report.pdf
  lms_r = rgb2lms(1,:);
  lms_b = rgb2lms(3,:);
  
  % Here we take the white points of it. We just multiply it by ones
  lms_w = rgb2lms * ones(3, 1);

  % After this is done, we need to get the cross section between the spaces we just got
  %   Again, this is explained in the PDF file
  p0 = cross(lms_w, lms_b');
  p1 = cross(lms_w, lms_r');

  % These are the transform matrices to perform the conversion to:
  %   Protanopia -> L cone defective (Long wave-length cones)
  %   Deuteranopia -> M cone defective (Medium wave-length cones)
  %   Tritanopia -> S cone defective (Small wave-length cones)
  %
  % The calculation is provided in the Android SDK
  %   https://android.googlesource.com/platform/frameworks/native/+/14e8b01/services/surfaceflinger/Effects/Daltonizer.cpp
  lms2lms_p = [  0.0000       0.0000,  0.0000;
               -p0(2)/p0(1),  1.0000,  0.0000;
               -p0(3)/p0(1),  0.0000,  1.0000]';
                    
  lms2lms_d = [ 1.0000, -p0(1)/p0(2),  0.0000;
                0.0000,       0.0000,  0.0000;
                0.0000, -p0(3)/p0(2),  1.0000]';

  lms2lms_t = [ 1.0000, 0.0000, -p1(1)/p1(3);
                0.0000, 1.0000, -p1(2)/p1(3);
                0.0000, 0.0000,      0.0000]';

  % Now we just apply the requested transformations on every pixel
  for i = 1:sizeRGB(1)
      for j = 1:sizeRGB(2)
          % First, we transform the RGB space to LMS space. This transformation varies depending on
          %   the transformation previously chosen.
          rgb = RGB(i,j,:);
          rgb = rgb(:);
          
          LMS(i,j,:) = rgb2lms * rgb;
          lms = LMS(i,j,:);
          lms = lms(:);
          
          % Then, we transform the current pixel to its correspondence within the colorblind system
          % So, in this step, we are transforming the pixels to colorblind LMS values
          LMSp(i,j,:) = lms2lms_p * lms;
          LMSd(i,j,:) = lms2lms_d * lms;
          LMSt(i,j,:) = lms2lms_t * lms;
          
          % Finally, we transform back to RGB space in order to show it
          lmsp = LMSp(i,j,:);
          lmsp = lmsp(:);
          
          lmsd = LMSd(i,j,:);
          lmsd = lmsd(:);
          
          lmst = LMSt(i,j,:);
          lmst = lmst(:);
          
          RGBp(i,j,:) = lms2rgb * lmsp;
          RGBd(i,j,:) = lms2rgb * lmsd;
          RGBt(i,j,:) = lms2rgb * lmst;
      end
  end

  % Before we finish the process we normalize the pixels value
  RGBp = uint8(RGBp);
  RGBd = uint8(RGBd);
  RGBt = uint8(RGBt);

  % Finally we write the images into separate files
  imwrite(RGBp,['_p.jpeg'],'jpeg');
  imwrite(RGBd,['_d.jpeg'],'jpeg');
  imwrite(RGBt,['_t.jpeg'],'jpeg');
  imwrite(BW, ['_bw.jpeg'], 'jpeg');
endfunction

%http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
function rgb2xyz = get_rgb2xyz_matrix(type)
  switch type
    case "adobe_rgb_d65"
      rgb2xyz = [ 0.5767309  0.1855540  0.1881852;
                  0.2973769  0.6273491  0.0752741;
                  0.0270343  0.0706872  0.9911085];
    case "adobe_rgb_d50"
      rgb2xyz = [ 0.6097559  0.2052401  0.1492240;
                  0.3111242  0.6256560  0.0632197;
                  0.0194811  0.0608902  0.7448387];
    case "apple_rgb_d65"
      rgb2xyz = [ 0.4497288  0.3162486  0.1844926;
                  0.2446525  0.6720283  0.0833192;
                  0.0251848  0.1411824  0.9224628];
    case "apple_rgb_d50"
      rgb2xyz = [ 0.4755678  0.3396722  0.1489800;
                  0.2551812  0.6725693  0.0722496;
                  0.0184697  0.1133771  0.6933632];
    case "best_rgb"
      rgb2xyz = [ 0.6326696  0.2045558  0.1269946;
                  0.2284569  0.7373523  0.0341908;
                  0.0000000  0.0095142  0.8156958];
    case "beta_rgb"
      rgb2xyz = [ 0.6712537  0.1745834  0.1183829;
                  0.3032726  0.6637861  0.0329413;
                  0.0000000  0.0407010  0.7845090];
    case "bruce_rgb_d65"
      rgb2xyz = [ 0.4674162  0.2944512  0.1886026;
                  0.2410115  0.6835475  0.0754410;
                  0.0219101  0.0736128  0.9933071];
    case "bruce_rgb_d50"
      rgb2xyz = [ 0.4941816  0.3204834  0.1495550;
                  0.2521531  0.6844869  0.0633600;
                  0.0157886  0.0629304  0.7464909];
    case "cie_rgb_e"
      rgb2xyz = [ 0.4887180  0.3106803  0.2006017;
                  0.1762044  0.8129847  0.0108109;
                  0.0000000  0.0102048  0.9897952];
    case "cie_rgb_d50"
      rgb2xyz = [ 0.6343706  0.1852204  0.1446290;
                  0.3109496  0.5915984  0.0974520;
                 -0.0011817  0.0555518  0.7708399];
    case "color_match_rgb"
      rgb2xyz = [ 0.5093439  0.3209071  0.1339691;
                  0.2748840  0.6581315  0.0669845;
                  0.0242545  0.1087821  0.6921735];
    case "don_rgb_4"
      rgb2xyz = [ 0.6457711  0.1933511  0.1250978;
                  0.2783496  0.6879702  0.0336802;
                  0.0037113  0.0179861  0.8035125];
    case "eci_rgb"
      rgb2xyz = [ 0.6502043  0.1780774  0.1359384;
                  0.3202499  0.6020711  0.0776791;
                  0.0000000  0.0678390  0.7573710];
    case "ekta_space_ps5"
      rgb2xyz = [ 0.5938914  0.2729801  0.0973485;
                  0.2606286  0.7349465  0.0044249;
                  0.0000000  0.0419969  0.7832131];
    case "ntsc_rgb_c"
      rgb2xyz = [ 0.6068909  0.1735011  0.2003480;
                  0.2989164  0.5865990  0.1144845;
                  0.0000000  0.0660957  1.1162243];
    case "ntsc_rgb_d50"
      rgb2xyz = [ 0.6343706  0.1852204  0.1446290;
                  0.3109496  0.5915984  0.0974520;
                 -0.0011817  0.0555518  0.7708399];
    case "pal_secam_rgb_d65"
      rgb2xyz = [ 0.4306190  0.3415419  0.1783091;
                  0.2220379  0.7066384  0.0713236;
                  0.0201853  0.1295504  0.9390944];
    case "pal_secam_rgb_d50"
      rgb2xyz = [ 0.4552773  0.3675500  0.1413926;
                  0.2323025  0.7077956  0.0599019;
                  0.0145457  0.1049154  0.7057489];
    case "pro_photo_rgb"
      rgb2xyz = [ 0.7976749  0.1351917  0.0313534;
                  0.2880402  0.7118741  0.0000857;
                  0.0000000  0.0000000  0.8252100];
    case "smptec_rgb_d65"
      rgb2xyz = [ 0.3935891  0.3652497  0.1916313;
                  0.2124132  0.7010437  0.0865432;
                  0.0187423  0.1119313  0.9581563];
    case "smptec_rgb_d50"
      rgb2xyz = [ 0.4163290  0.3931464  0.1547446;
                  0.2216999  0.7032549  0.0750452;
                  0.0136576  0.0913604  0.7201920];
    case "srgb_d65"
      rgb2xyz = [ 0.4124564  0.3575761  0.1804375;
                  0.2126729  0.7151522  0.0721750;
                  0.0193339  0.1191920  0.9503041];
    case "srgb_d50"
      rgb2xyz = [ 0.4360747  0.3850649  0.1430804;
                  0.2225045  0.7168786  0.0606169;
                  0.0139322  0.0971045  0.7141733];
    case "wide_gamut_rgb"
      rgb2xyz = [	0.7161046  0.1009296  0.1471858;
                  0.2581874  0.7249378  0.0168748;
                  0.0000000  0.0517813  0.7734287];
    otherwise
      rgb2xyz = 0
  endswitch
endfunction

%https://en.wikipedia.org/wiki/LMS_color_space#XYZ_to_LMS
function xyz2lms = get_xyz2lms_matrix(type)
  switch type
    case "equal_energy"
      xyz2lms = [ 0.38971  0.68898 -0.07868; 
                 -0.22981  1.18340  0.04641;
                  0.00000  0.00000  1.00000];
    case "equal_energy_d65"
      xyz2lms = [ 0.4002  0.7076 -0.0808;
                 -0.2263  1.1653  0.0457;
                  0.0000  0.0000  0.9182];
    case "ciecam97"
      xyz2lms = [ 0.8951  0.2664 -0.1614;
                 -0.7502  1.7135  0.0367;
                  0.0389 -0.0685  1.0296];
    case "ciecam02"
      xyz2lms = [ 0.7328 0.4296 -0.1624;
                 -0.7036 1.6975  0.0061;
                  0.0030 0.0136  0.9834];
    otherwise
      xyz2lms = 0
  endswitch
endfunction