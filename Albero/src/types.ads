with Ada.Text_IO; use Ada.Text_IO;
with Ada.Float_Text_IO; use Ada.Float_Text_IO;
with Ada.Numerics; use Ada.Numerics;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;

package types is
   type Longitudine is new Float; -- range -180.0 .. 180.0;
   type Latitudine is new Float; -- range -90.0 .. 90.0;
   type altitudine is new Float; -- range -180.0 .. 180.0;
   subtype Theta is Float; -- range 0.0 .. 360.0;
   
   type Point is record
      long : Longitudine;
      lat : Latitudine;
      alt : altitudine := 5000.0;
   end record;
      
end types;
