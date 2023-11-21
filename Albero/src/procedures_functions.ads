with Ada.Text_IO; use Ada.Text_IO;
with Ada.Float_Text_IO; use Ada.Float_Text_IO;
with Ada.Numerics; use Ada.Numerics;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;
with types; use types;

package Procedures_functions is
   function Deg2Rad (pippo:in out float) return float;
   function Rad2Deg (pippo: in out float) return float;
   procedure CalculatePoint(A:in Point; Teta:in Theta; dist:in Float; B:out Point);
   procedure CalculatePointHaversine(A:in Point; Teta:in Theta; dist:in Float; B:out Point);
   procedure CalculatePointVincenty(A:in Point; Teta:in Theta; dist:in Float; B:out Point);
   procedure DisplayPoint(P:in Point; Name:in String);
   function GetUserInput(pippo : String) return Float;
   function GetUserInputInRange(pippo:String; Min, Max:Float) return Float;
   function GetUserLatitudeInput(pippo:String) return Latitudine;
   function GetUserLongitudeInput(pippo:String) return Longitudine;
   function GetUserThetaInput(pippo:String) return Theta;
   function ViewPoint(P:in Point) return String;
   function remove_spaces(actual: String) return String;

end Procedures_functions;
