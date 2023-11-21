with Ada.Text_IO; use Ada.Text_IO;
with Ada.Float_Text_IO; use Ada.Float_Text_IO;
with Ada.Integer_Text_IO; use Ada.Integer_Text_IO;
with Ada.Numerics; use Ada.Numerics;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;
with Procedures_functions; use Procedures_functions;
with Ada.Strings.Unbounded; use Ada.Strings.Unbounded;
with types; use types;

procedure Main is
   Origin : Point;
   Teta : Theta;
   dist : Float;
   B, C, D, E, F, G, H, I, L, M, N : Point;
   Num : Integer :=0;

   Ft : File_Type;
   Fg: File_Type;
   File_Name: Unbounded_String := Null_Unbounded_String;
   piatta_file : String:= "terrapiatta.kml";
   sferica_file: String:= "sferica.kml";
   vinc_file: String := "vincenty.kml";
   General_File: constant String:= "general.kml";
   wanted_len : Integer;
   wanted_str : String := "<coordinates>";
   str_bed: String (1..13);
   trial_bed: String(1..18);
   tmp_unbounded_str: Unbounded_String := Null_Unbounded_String;
   tmp_char: String := "A";
   spaced_str : Unbounded_String := Null_Unbounded_String;
   all_info : Unbounded_String := Null_Unbounded_String;

begin
   -- Input user for the origine point
   Origin.long := GetUserLongitudeInput("Insert the longitude of the first point (A) < -180 .. 180 >: ");
   Origin.lat := GetUserLatitudeInput("Insert the latitude of the first point (A) < -90 .. 90 >: ");

   -- Input user for the angle Theta
   Teta := GetUserThetaInput("Insert the angle Theta in grad < 0 .. 360 >: ");
   Teta := Deg2Rad(Teta);

   -- Input user for the distance
   dist := GetUserInput("Insert the distance in miglia: ");
   New_Line;

   while (Num not in 1..3) loop
      Put_Line("Insert 1 for the 2D Earth Model");
      Put_Line("Insert 2 for the Spheric Earth Model");
      Put_Line("Insert 3 for the WGS84 Earth Model");
      Get (Num);

      case Num is
         when 1 =>
            -- Calcolo dei punti
            CalculatePoint(Origin, Teta + (Pi/2.0), dist/60.0, B);
            CalculatePoint(B, Teta, dist/60.0, C);
            CalculatePoint(C, Teta + (Pi/2.0), dist/60.0, D);
            CalculatePoint(D, Teta - (Pi/4.0), dist/60.0 * Sqrt(2.0), E);
            CalculatePoint(E, Teta + (Pi/2.0), dist/60.0, F);
            CalculatePoint(F, Teta - (Pi/4.0), (3.0 / 2.0) * dist/60.0 * Sqrt(2.0), G);
            CalculatePoint(G, Teta - (3.0*Pi/4.0), (3.0 / 2.0) * dist/60.0 * Sqrt(2.0), H);
            CalculatePoint(H, Teta + (Pi/2.0), dist/60.0, I);
            CalculatePoint(I, Teta - (3.0*Pi/4.0), dist/60.0 * Sqrt(2.0), L);
            CalculatePoint(L, Teta + (Pi/2.0), dist/60.0, M);
            CalculatePoint(M, Teta + Pi, dist/60.0, N);
            File_Name := To_Unbounded_String(piatta_file);

            -- N := Origin;
         when 2 =>
            -- Calcolo dei punti
            CalculatePointHaversine(Origin, Teta + (Pi/2.0), dist, B);
            CalculatePointHaversine(B, Teta, dist, C);
            CalculatePointHaversine(C, Teta + (Pi/2.0), dist, D);
            CalculatePointHaversine(D, Teta - (Pi/4.0), dist * Sqrt(2.0), E);
            CalculatePointHaversine(E, Teta + (Pi/2.0), dist, F);
            CalculatePointHaversine(F, Teta - (Pi/4.0), (3.0 / 2.0) * dist * Sqrt(2.0), G);
            CalculatePointHaversine(G, Teta - (3.0*Pi/4.0), (3.0 / 2.0) * dist * Sqrt(2.0), H);
            CalculatePointHaversine(H, Teta + (Pi/2.0), dist, I);
            CalculatePointHaversine(I, Teta - (3.0*Pi/4.0), dist * Sqrt(2.0), L);
            CalculatePointHaversine(L, Teta + (Pi/2.0), dist, M);
            CalculatePointHaversine(M, Teta + Pi, dist, N);
            File_Name := To_Unbounded_String(sferica_file);

         when 3 =>
            -- Calcolo dei punti
            CalculatePointVincenty(Origin, Teta + (Pi/2.0) , dist, B);
            CalculatePointVincenty(B, Teta, dist, C);
            CalculatePointVincenty(C, Teta + (Pi/2.0) , dist, D);
            CalculatePointVincenty(D, Teta - (Pi/4.0) , dist * Sqrt(2.0), E);
            CalculatePointVincenty(E, Teta + (Pi/2.0) , dist, F);
            CalculatePointVincenty(F, Teta - (Pi/4.0) , (3.0 / 2.0) * dist * Sqrt(2.0), G);
            CalculatePointVincenty(G, Teta - (3.0*Pi/4.0) , (3.0 / 2.0) * dist * Sqrt(2.0), H);
            CalculatePointVincenty(H, Teta + (Pi/2.0) , dist, I);
            CalculatePointVincenty(I, Teta  - (3.0*Pi/4.0) , dist * Sqrt(2.0), L);
            CalculatePointVincenty(L, Teta + (Pi/2.0) , dist, M);
            CalculatePointVincenty(M, Teta + Pi , dist, N);
            File_Name := To_Unbounded_String(vinc_file);

         when others =>
            Put_Line("You have to insert a valid integer from 1 to 3");

      end case;
   end loop;


   -- Visualizzazione dei punti
   DisplayPoint(Origin, "A");
   DisplayPoint(B, "B");
   DisplayPoint(C, "C");
   DisplayPoint(D, "D");
   DisplayPoint(E, "E");
   DisplayPoint(F, "F");
   DisplayPoint(G, "G");
   DisplayPoint(H, "H");
   DisplayPoint(I, "I");
   DisplayPoint(L, "L");
   DisplayPoint(M, "M");
   DisplayPoint(N, "N");

   Create(Ft, Out_File, "C:\Users\ibrahim\Desktop\Albero\" & To_String(File_Name));
   --Reset(Ft);   -- Reset file you are going to write in
   Open(Fg, In_File, "C:\Users\ibrahim\Desktop\Albero\" & General_File); -- Read from general file
   while not End_Of_File(Fg) loop
      tmp_unbounded_str := To_Unbounded_String(Get_Line(Fg));
      Put_Line(Ft, To_String(tmp_unbounded_str));
      wanted_len := Length(tmp_unbounded_str);
      if wanted_len = 18 then
         trial_bed := To_String(tmp_unbounded_str);
         for I in 1..13 loop
            str_bed(I) := trial_bed(I+5);
         end loop;
         if str_bed = wanted_str then
            Put_Line(Ft, remove_spaces(ViewPoint(Origin)));
            Put_Line(Ft, remove_spaces(ViewPoint(B)));
            Put_Line(Ft, remove_spaces(ViewPoint(C)));
            Put_Line(Ft, remove_spaces(ViewPoint(D)));
            Put_Line(Ft, remove_spaces(ViewPoint(E)));
            Put_Line(Ft, remove_spaces(ViewPoint(F)));
            Put_Line(Ft, remove_spaces(ViewPoint(G)));
            Put_Line(Ft, remove_spaces(ViewPoint(H)));
            Put_Line(Ft, remove_spaces(ViewPoint(I)));
            Put_Line(Ft, remove_spaces(ViewPoint(L)));
            Put_Line(Ft, remove_spaces(ViewPoint(M)));
            Put_Line(Ft, remove_spaces(ViewPoint(N)));
         end if;
      end if;
   end loop;
   Close(Fg);
   Close(Ft);
end Main;
