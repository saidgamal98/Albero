with Ada.Text_IO; use Ada.Text_IO;
with Ada.Float_Text_IO; use Ada.Float_Text_IO;
with Ada.Numerics; use Ada.Numerics;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;
with types; use types;
with ada.strings.Unbounded; use ada.Strings.Unbounded;

package body Procedures_functions is
      
   -- Function to convert from grads to rads
   function Deg2Rad (pippo: in out Float) return Float is
   begin
      pippo := pippo*Pi/180.0;
        return pippo;
   end Deg2Rad;
   
   function Rad2Deg (pippo: in out Float) return Float is
   begin
      pippo := pippo*180.0/Pi;
        return pippo;
   end Rad2Deg;
   
   -- Procedure that calculate the earth 2D model
   procedure CalculatePoint(A:in Point; Teta:in Theta; dist:in Float; B:out Point) is
   begin
      B.long := A.long + Longitudine(dist * Sin(Float(Teta)));
      if B.long > Longitudine(180) then
         B.long := B.long - Longitudine(360);
      elsif B.long < Longitudine(-180) then
         B.long := B.long + Longitudine(360);
      end if;   
      B.lat := A.lat + Latitudine(dist * Cos(Float(Teta)));
      if B.lat > Latitudine(90) then
         B.lat :=  Latitudine(180) - B.lat;
      elsif B.lat < Latitudine(-90) then
         B.lat := Latitudine(-180) - B.lat;
      end if;   
   end CalculatePoint;
   
   -- Procedure that calculate the earth sphere model
   procedure CalculatePointHaversine(A:in Point; Teta:in Theta; dist:in Float; B:out Point) is
      d : Float := dist/3444.0;
      l1 : Float;
      temp_lat : float; 
      temp_lon : float; 
   begin
      l1 := Float(A.lat);
      l1 := Deg2Rad(l1);
      temp_lat := Arcsin(Sin(l1) * Cos(d) + 
                       Cos(l1) * Sin(d) * Cos(Teta));
      B.lat := Latitudine(Rad2Deg(temp_lat));
      -- Control on the range of the latitude
      if B.lat > Latitudine(90) then
         B.lat :=  Latitudine(180) - B.lat;
      elsif B.lat < Latitudine(-90) then
         B.lat := Latitudine(-180) - B.lat;
      end if; 
  
         temp_lon := float(A.long); 
         temp_lon := Deg2Rad(temp_lon) + Arctan((Cos(l1) * Sin(d) * Sin(Float(Teta))), 
                                        Cos(d) - (Sin(l1) * Sin(temp_lat)));
         B.long := Longitudine(Rad2Deg(temp_lon)); 
      
      -- Control on the range of the latitude
      if B.long > Longitudine(180) then
         B.long := B.long - Longitudine(360);
      elsif B.long < Longitudine(-180) then
         B.long := B.long + Longitudine(360);
      end if;     
   end CalculatePointHaversine;
   
   -- Procedure to calculate the destination point given start point latitude / longitude (numeric degrees),
   -- bearing (numeric degrees) and distance (in Nautical miles).
   procedure CalculatePointVincenty(A:in Point; Teta:in Theta; dist:in Float; B:out Point) is
      a1 : Float := 3440.065; --6378137.0; --m
      b1 : Float := 3432.505;-- 6356752.3142; -- m
      f : Float ; -- 1.0 / 298.257;
      s : Float := dist;--*1852.0; --m
      alpha1 : Float := Teta;
      sin_alpha1 : Float;
      cos_alpha1 : Float;
      A_lat_Deg : Float;
      A_long_Deg : Float;
      A_lat_toRad : Float;
      A_long_toRad : Float;
      tan_U1 : Float;
      cos_U1 : Float;
      sin_U1 : Float;
      sigma1 : Float;
      sin_alpha : Float;
      cos_Sq_alpa : Float;
      u_Sq : Float;
      AA : Float;
      BB : Float;
      sigma : Float;
      sigmaP : Float;
      cos2sigmaM : Float;
      sin_sigma : Float;
      cos_sigma : Float;
      delta_sigma : Float;
      temp : Float;
      B_lat_Rad : Float;
      B_lat_toDeg : Float;
      lambda : Float;
      CC : Float;
      B_long_Rad : Float;
      B_long_toDeg : Float;
      DEBUG : Boolean := True;
      
   begin
      f := (a1 - b1) / a1; --1.0 / 298.257;   --2.19765E-03;
      sin_alpha1 := Sin(alpha1);
      cos_alpha1 := Cos(alpha1);
      A_lat_Deg := Float(A.lat);
      A_lat_toRad := Deg2Rad(A_lat_Deg);
      A_long_Deg := Float(A.long);
      A_long_toRad := Deg2Rad(A_long_Deg);
      tan_U1 := (1.0 - f) * Tan(A_lat_toRad);
      cos_U1 := 1.0 / Sqrt(1.0 + tan_U1 * tan_U1);
      sin_U1 := tan_U1 * cos_U1;
      sigma1 := Arctan(tan_U1/cos_alpha1);
      sin_alpha := cos_U1 * sin_alpha1;
      cos_Sq_alpa := 1.0 - sin_alpha * sin_alpha;
      u_Sq := cos_Sq_alpa * (a1*a1 - b1*b1) / (b1*b1);
      AA := 1.0 + u_Sq / 16384.0 * (4096.0 + u_Sq * (-768.0 + u_Sq * (320.0 - 175.0 * u_Sq)));
      BB := u_Sq / 1024.0 * (256.0 + u_Sq * (-128.0 + u_Sq * (74.0 - 47.0 * u_Sq)));
      sigma := s / (b1 * AA);
      sigmaP := 0.0;
      
      while (abs(sigma - sigmaP) > 1.0E-12) loop
         cos2sigmaM := Cos(2.0*sigma1 + sigma);
         sin_sigma := Sin(sigma);
         cos_sigma := Cos(sigma);
         delta_sigma := BB * sin_sigma * (cos2sigmaM + BB / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos2sigmaM 
                                          * cos2sigmaM) - BB / 6.0 * cos2sigmaM * (-3.0 + 4.0 * sin_sigma 
                                            * sin_sigma) * (-3.0 + 4.0 * cos2sigmaM * cos2sigmaM)));
         sigmaP := sigma;
         sigma := s / (b1 * AA) + delta_sigma;
      end loop;
      
      temp := sin_U1 * sin_sigma - cos_U1 * cos_sigma * cos_alpha1;
      B_lat_Rad := Arctan(sin_U1 * cos_sigma + cos_U1 * sin_sigma * cos_alpha1/((1.0 - f) * Sqrt(sin_alpha * sin_alpha + temp * temp)));
      B_lat_toDeg := Rad2Deg(B_lat_Rad);
      B.lat := Latitudine(B_lat_toDeg);
      lambda := Arctan(sin_sigma * sin_alpha1/(cos_U1 * cos_sigma - sin_U1 * sin_sigma * cos_alpha1));
      CC := f / 16.0 * cos_Sq_alpa * (4.0 + f * (4.0 - 3.0 * cos_Sq_alpa));
      B_long_Rad := lambda - (1.0 - CC) * f * sin_alpha * (sigma + CC * sin_sigma * (cos2sigmaM + CC * cos_sigma * (-1.0 + 2.0 * cos2sigmaM * cos2sigmaM)));
      --B_long_toDeg := Rad2Deg(A_lat_toRad + B_long_Rad);
      --B.long := Longitudine(B_long_toDeg);
      B_long_toDeg := Rad2Deg(B_long_Rad);
      B.long := A.long + Longitudine(B_long_toDeg);
      -- Teta := Arctan(sin_alpha, -temp); -- final bearing
      
      if DEBUG then
            Put_Line("Flattening: " & Float'Image(f));
            Put_Line("Sin of bearing: " & Float'Image(sin_alpha1));
            Put_Line("Cos of bearing: " & Float'Image(cos_alpha1));
            Put_Line("Tan U: " & Float'Image(tan_U1));
            Put_Line("Cos U: " & Float'Image(cos_U1));
            Put_Line("Sin U: " & Float'Image(sin_U1));
            Put_Line("Ang_dist: " & Float'Image(sigma1));
            Put_Line("Sin_a: " & Float'Image(sin_alpha));
            Put_Line("Cos2a: " & Float'Image(cos_Sq_alpa));
            Put_Line("U2: " & Float'Image(u_Sq));
            Put_Line("A: " & Float'Image(AA));
            Put_Line("B: " & Float'Image(BB));
            Put_Line("Absolute: " & Float'Image(abs(sigma-sigmaP)));
            --Put_Line("New Lat : " & Float'Image());
            --Put_Line("New Longitude " & Float'Image());
         end if;
   end CalculatePointVincenty;
   
   
   procedure DisplayPoint(P:in Point; Name:in String) is
   begin
      -- Put_Line(Name & ": (" & P.long'Image & ", " & P.lat'Image & ", " & P.alt'Image & ")");
      Put_Line(P.long'Image & "," & P.lat'Image & "," & P.alt'Image);
   end DisplayPoint;
   
   function ViewPoint(P:in Point) return String is
   begin
      return P.long'Image & "," & P.lat'Image & "," & P.alt'Image;
   end ViewPoint;
   
   -- Function to get a value from the user 
   function GetUserInput(pippo:String) return Float is
      Input_Value : Float;
   begin
      Put(pippo);
      Get(Input_Value);
      return Input_Value;
   end GetUserInput;
   
   -- Function to get a value from the user in a certain range
   function GetUserInputInRange(pippo:String; Min, Max:Float) return Float is
      User_Input : Float;
   begin
      loop
         User_Input := GetUserInput(pippo);
         exit when User_Input >= Min and User_Input <= Max;
         Put_Line("Input not valid. Put a value between " & Min'Image & " and " & Max'Image);
      end loop;
      return User_Input;
   end GetUserInputInRange;
  
   -- Function to get a value from the user for the longitude
   function GetUserLongitudeInput(pippo:String) return Longitudine is
   begin
      return Longitudine(GetUserInputInRange(pippo, -180.0, 180.0));
   end GetUserLongitudeInput;
   
   -- Function to get a value from the user for the latitude
   function GetUserLatitudeInput(pippo:String) return Latitudine is
   begin
      return Latitudine(GetUserInputInRange(pippo, -90.0, 90.0));
   end GetUserLatitudeInput;

   -- Function to get a value from the user for the angle Theta
   function GetUserThetaInput(pippo:String) return Theta is
   begin
      return Theta(GetUserInputInRange(pippo, 0.0, 360.0));
   end GetUserThetaInput;
   
   function remove_spaces(actual: String) return String is
      no_space_str : Unbounded_String :=Null_Unbounded_String;
   begin
      for I in actual'range loop
         if actual(I) /= ' ' then
            Append(no_space_str, actual(I));
         end if;
      end loop;
      return To_String(no_space_str);
   end remove_spaces;
   
end Procedures_functions;
