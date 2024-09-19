function result = Rxx_gaussian(d_norm,phi_0_deg,sigma_deg,delta_phi_deg)
% result = Rxx_gaussian(d_norm,phi_0_deg,sigma_deg,delta_phi_deg)
%
% Computes the correlation of the real and imaginary parts of the
% signals in the case of a truncated gaussian Power Azimuth
% Spectrum (PAS) at spacings given by d_norm. The PAS is
% characterised by the Angle Of Arrival (AOA) phi_0_deg and by its
% standard deviation sigma_deg.
%
% This program calls a procedure named erfcomp. It computes the
% error function of a complex argument according to Abramowitz,
% "Error Functions and Fresnel Integrals", p. 299.
%
%
% STANDARD DISCLAIMER
%
% CSys is furnishing this item "as is". CSys does not provide any
% warranty of the item whatsoever, whether express, implied, or
% statutory, including, but not limited to, any warranty of
% merchantability or fitness for a particular purpose or any
% warranty that the contents of the item will be error-free.
%
% In no respect shall CSys incur any liability for any damages,
% including, but limited to, direct, indirect, special, or
% consequential damages arising out of, resulting from, or any way
% connected to the use of the item, whether or not based upon
% warranty, contract, tort, or otherwise; whether or not injury was
% sustained by persons or property or otherwise; and whether or not
% loss was sustained from, or arose out of, the results of, the
% item, or any services that may be provided by CSys.
%
% (c) Laurent Schumacher, AAU-TKN/IES/KOM/CPK/CSys - July 2001

D               = 2*pi*d_norm;
% Conversion degree -> rad
phi_0_rad       = phi_0_deg*(pi/180);
sigma_rad       = sigma_deg*(pi/180);
delta_phi_rad   = delta_phi_deg*(pi/180);

m               = 0;
result          = zeros(size(D));
tmp_xx_gaussian = ones(size(D));
while (m < 100)
    m = m +1;
    A = erfcomp((delta_phi_rad/(sigma_rad*sqrt(2)))-(j*m*sigma_rad*sqrt(2)))...
      - erfcomp(((-1)*delta_phi_rad/(sigma_rad*sqrt(2)))-(j*m*sigma_rad*sqrt(2)));
    if (isnan(A))
        disp('Warning: NaN reached in Rxx_gaussian - Breaking loop');
        break;
    end;
    tmp_xx_gaussian = besselj(2*m,D).*cos(2*m*phi_0_rad).*exp((-2)* ...
						  (sigma_rad^2)*(m^2)).*real(A);
    result          = result + tmp_xx_gaussian;
end;