using System;
using System.IO;
using System.Text.Json;
using MeanlineSolver.Core;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("Reading design.json...\n");

        string json = File.ReadAllText("design.json");
        MeanlineInputs inputs = JsonSerializer.Deserialize<MeanlineInputs>(json);

        if (inputs == null)
        {
            Console.WriteLine("ERROR: JSON parse edilemedi.");
            return;
        }

        var solver = new MeanlineSolverMain();
        var results = solver.Solve(inputs);

        foreach (var result in results)
        {
            PrintResult(result);
        }

        Console.WriteLine("\nFinished.");
    }

    static void PrintResult(StageResult r)
    {
        Console.WriteLine("\n=========== STAGE OUTPUT (PDF FORMAT) ===========");

        Console.WriteLine("\n--- Flow Inlet Angles ---");
        Console.WriteLine($"Alpha1 (flow)         : {RadToDeg(r.Station1.Alpha):F3} deg");
        Console.WriteLine($"Beta1 (relative)      : {RadToDeg(r.Station1.Beta):F3} deg");

        Console.WriteLine("\n--- Flow Outlet Angles ---");
        Console.WriteLine($"Alpha2 (flow)         : {RadToDeg(r.Station2.Alpha):F3} deg");
        Console.WriteLine($"Beta2 (relative)      : {RadToDeg(r.Station2.Beta):F3} deg");
        Console.WriteLine($"Alpha3 (stator exit)  : {RadToDeg(r.Station3.Alpha):F3} deg");

        Console.WriteLine("\n--- Blade Angles ---");
        Console.WriteLine($"Rotor Inlet Blade     : {RadToDeg(r.RotorInletAngle):F3} deg");
        Console.WriteLine($"Rotor Exit Blade      : {RadToDeg(r.RotorExitAngle):F3} deg");
        Console.WriteLine($"Stator Inlet Blade    : {RadToDeg(r.StatorInletAngle):F3} deg");
        Console.WriteLine($"Stator Exit Blade     : {RadToDeg(r.StatorExitAngle):F3} deg");

        Console.WriteLine("\n--- Velocities ---");
        Console.WriteLine($"Ca1, Ct1, C1          : {r.Station1.Ca}, {r.Station1.Ct}, {r.Station1.C}");
        Console.WriteLine($"Ca2, Ct2, C2          : {r.Station2.Ca}, {r.Station2.Ct}, {r.Station2.C}");
        Console.WriteLine($"Ca3, Ct3, C3          : {r.Station3.Ca}, {r.Station3.Ct}, {r.Station3.C}");

        Console.WriteLine("\n--- Incidence & Deviation ---");
        Console.WriteLine($"Incidence             : {RadToDeg(r.Incidence):F3} deg");
        Console.WriteLine($"Deviation             : {RadToDeg(r.Deviation):F3} deg");

        Console.WriteLine("\n--- Camber Angles ---");
        Console.WriteLine($"Rotor Camber          : {RadToDeg(r.CamberAngleRotor):F3} deg");
        Console.WriteLine($"Stator Camber         : {RadToDeg(r.CamberAngleStator):F3} deg");

        Console.WriteLine("\n--- Chord Lengths ---");
        Console.WriteLine($"Rotor Chord           : {r.ChordLengthRotor}");
        Console.WriteLine($"Stator Chord          : {r.ChordLengthStator}");

        Console.WriteLine("\n--- Solidities ---");
        Console.WriteLine($"Solidity (final)      : {r.Solidity}");

        Console.WriteLine("\n--- Radiuses ---");
        Console.WriteLine($"Station1 (rh, rm, rt) : {r.Station1.HubRadius:F4}, {r.Station1.Radius:F4}, {r.Station1.TipRadius:F4}");
        Console.WriteLine($"Station2 (rh, rm, rt) : {r.Station2.HubRadius:F4}, {r.Station2.Radius:F4}, {r.Station2.TipRadius:F4}");
        Console.WriteLine($"Station3 (rh, rm, rt) : {r.Station3.HubRadius:F4}, {r.Station3.Radius:F4}, {r.Station3.TipRadius:F4}");

        Console.WriteLine("\n--- Flow Properties ---");
        Console.WriteLine($"T1, P1, rho1          : {r.Station1.T}, {r.Station1.P}, {r.Station1.Rho}");
        Console.WriteLine($"T2, P2, rho2          : {r.Station2.T}, {r.Station2.P}, {r.Station2.Rho}");
        Console.WriteLine($"T3, P3, rho3          : {r.Station3.T}, {r.Station3.P}, {r.Station3.Rho}");

        Console.WriteLine("\n--- Performance ---");
        Console.WriteLine($"Δh0                   : {r.DeltaH0}");
        Console.WriteLine($"Pressure Ratio        : {r.PressureRatio}");
        Console.WriteLine($"Efficiency            : {r.Efficiency}");

        Console.WriteLine("==================================================\n");
    }

    static double RadToDeg(double rad) => rad * 180.0 / Math.PI;
}
