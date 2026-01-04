using System;
using System.IO;
using System.Text.Json;
using MeanlineSolver.Core;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("Reading design.json...\n");

        // 1) JSON dosyasını oku
        string json = File.ReadAllText("design.json");

        // 2) JSON → C# modele dönüştür
        MeanlineInputs inputs = JsonSerializer.Deserialize<MeanlineInputs>(json);

        if (inputs == null)
        {
            Console.WriteLine("ERROR: JSON parse edilemedi.");
            return;
        }

        if (inputs.Stages == null || inputs.Stages.Count == 0)
        {
            Console.WriteLine("ERROR: JSON içinde hiç 'Stages' bulunamadı.");
            return;
        }

        // 3) Stage solver oluştur
        var solver = new MeanlineSolverMain();
        var results = solver.Solve(inputs);

        foreach (var result in results)
        {
            PrintResult(result);
        }


        Console.WriteLine("\nFinished.");
    }

    // 5) Çıktıları ekrana yazdıran fonksiyon
    static void PrintResult(StageResult r)
    {
        Console.WriteLine("\n=========== STAGE OUTPUT (PDF FORMAT) ===========");

        Console.WriteLine("\n--- Flow Inlet Angles ---");
        Console.WriteLine($"Alpha1 (flow)         : {r.Station1.Alpha}");
        Console.WriteLine($"Beta1 (relative)      : {r.Station1.Beta}");

        Console.WriteLine("\n--- Flow Outlet Angles ---");
        Console.WriteLine($"Alpha2 (flow)         : {r.Station2.Alpha}");
        Console.WriteLine($"Beta2 (relative)      : {r.Station2.Beta}");
        Console.WriteLine($"Alpha3 (stator exit)  : {r.Station3.Alpha}");

        Console.WriteLine("\n--- Blade Angles ---");
        Console.WriteLine($"Rotor Inlet Blade     : {r.RotorInletAngle}");
        Console.WriteLine($"Rotor Exit Blade      : {r.RotorExitAngle}");
        Console.WriteLine($"Stator Inlet Blade    : {r.StatorInletAngle}");
        Console.WriteLine($"Stator Exit Blade     : {r.StatorExitAngle}");

        Console.WriteLine("\n--- Velocities ---");
        Console.WriteLine($"Ca1, Ct1, C1          : {r.Station1.Ca}, {r.Station1.Ct}, {r.Station1.C}");
        Console.WriteLine($"Ca2, Ct2, C2          : {r.Station2.Ca}, {r.Station2.Ct}, {r.Station2.C}");
        Console.WriteLine($"Ca3, Ct3, C3          : {r.Station3.Ca}, {r.Station3.Ct}, {r.Station3.C}");

        Console.WriteLine("\n--- Incidence & Deviation ---");
        Console.WriteLine($"Incidence             : {r.Incidence}");
        Console.WriteLine($"Deviation             : {r.Deviation}");

        Console.WriteLine("\n--- Camber Angles ---");
        Console.WriteLine($"Rotor Camber          : {r.CamberAngleRotor}");
        Console.WriteLine($"Stator Camber         : {r.CamberAngleStator}");

        Console.WriteLine("\n--- Chord Lengths ---");
        Console.WriteLine($"Rotor Chord           : {r.ChordLengthRotor}");
        Console.WriteLine($"Stator Chord          : {r.ChordLengthStator}");

        Console.WriteLine("\n--- Solidities ---");
        Console.WriteLine($"Solidity (final)      : {r.Solidity}");

        Console.WriteLine("\n--- Radiuses ---");
        Console.WriteLine($"Radius Station1       : {r.Station1.Radius}");
        Console.WriteLine($"Radius Station2       : {r.Station2.Radius}");
        Console.WriteLine($"Radius Station3       : {r.Station3.Radius}");

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

}
