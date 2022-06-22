/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

#include "ImguiGPUControl.h"
#include "ImguiExtensions.h"
#include "imgui/imgui.h"
#include <imgui/imgui_internal.h>
#include <sstream>

#include "../../src/MolFlow.h"
extern MolFlow *mApp;

void ShowGPUWindow(MolFlow *mApp, bool *show_global_settings, bool &nbProcChanged, bool &recalcOutg,
                        bool &changeDesLimit, int &nbProc) {

    ImGui::PushStyleVar(
            ImGuiStyleVar_WindowMinSize,
            ImVec2(800.f, 0)); // Lift normal size constraint, however the presence of
    // a menu-bar will give us the minimum height we want.

    ImGui::Begin(
            "GPU Control", show_gpu,
            ImGuiWindowFlags_NoSavedSettings); // Pass a pointer to our bool
    // variable (the window will have
    // a closing button that will
    // clear the bool when clicked)
    ImGui::PopStyleVar(1);

    float gasMass = 2.4;
    if (ImGui::BeginTable("split", 2, ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableSetupColumn("Global settings");
        ImGui::TableSetupColumn("Simulation settings");
        ImGui::TableHeadersRow();
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::PushItemWidth(100);
        ImGui::InputDouble("Autosave frequency (minutes)",
                           &mApp->autoSaveFrequency, 0.00f, 0.0f, "%.1f");
        ImGui::PopItemWidth();
        ImGui::Checkbox(
                "Autosave only when simulation is running",
                reinterpret_cast<bool *>(
                        &mApp->autoSaveSimuOnly)); // Edit bools storing our window
        // open/close state
        ImGui::Checkbox(
                "Use .zip as default extension (otherwise .xml)",
                reinterpret_cast<bool *>(
                        &mApp->compressSavedFiles)); // Edit bools storing our window
        // open/close state
        ImGui::Checkbox(
                "Check for updates at startup",
                reinterpret_cast<bool *>(
                        &mApp->updateRequested)); // Edit bools storing our window
        // open/close state
        ImGui::Checkbox(
                "Auto refresh formulas",
                reinterpret_cast<bool *>(
                        &mApp->autoUpdateFormulas)); // Edit bools storing our window
        // open/close state
        ImGui::Checkbox("Anti-Aliasing",
                        &mApp->antiAliasing); // Edit bools storing our window
        // open/close state
        ImGui::Checkbox(
                "White Background",
                &mApp->whiteBg); // Edit bools storing our window open/close state
        ImGui::Checkbox("Left-handed coord. system",
                        &mApp->leftHandedView); // Edit bools storing our
        // window open/close state
        ImGui::Checkbox(
                "Highlight non-planar facets",
                &mApp->highlightNonplanarFacets); // Edit bools storing our window
        // open/close state
        ImGui::Checkbox("Highlight selected facets",
                        &mApp->highlightSelection); // Edit bools storing our
        // window open/close state
        ImGui::Checkbox("Use old XML format",
                        &mApp->useOldXMLFormat); // Edit bools storing our
        // window open/close state
        ImGui::TableNextColumn();
        ImGui::PushItemWidth(100);

        /* --- Simu settings ---*/
        static bool simChanged = false;
#if defined(MOLFLOW)
        static double gasMass = mApp->worker.model->wp.gasMass;
        static bool enableDecay = mApp->worker.model->wp.enableDecay;
        static double halfLife = mApp->worker.model->wp.halfLife;
        static bool lowFluxMode = mApp->worker.model->otfParams.lowFluxMode;
        static double lowFluxCutoff =
                mApp->worker.model->otfParams.lowFluxCutoff;

        simChanged =
                ImGui::InputRightSide("Gas molecular mass (g/mol)", &gasMass, "%g");
        simChanged = ImGui::Checkbox("", &enableDecay);
        if (!enableDecay) {
            ImGui::BeginDisabled();
        }
        ImGui::SameLine();
        simChanged = ImGui::InputRightSide("Gas half life (s)", &halfLife, "%g");

        if (!enableDecay) {
            ImGui::EndDisabled();
        }

        ImGui::BeginDisabled();

        // Use tmp var to multiply by 10
        double outgRate10 =
                mApp->worker.model->wp.finalOutgassingRate_Pa_m3_sec * 10.00;
        ImGui::InputRightSide("Final outgassing rate (mbar*l/sec)", &outgRate10,
                       "%.4g"); // 10: conversion Pa*m3/sec -> mbar*l/s
        ImGui::InputRightSide("Final outgassing rate (1/sec)",
                       &mApp->worker.model->wp.finalOutgassingRate,
                       "%.4g"); // In molecules/sec

        {
            char tmp[64];
            sprintf(tmp, "Tot.des. molecules [0 to %g s]",
                    mApp->worker.model->wp.latestMoment);
            ImGui::InputRightSide(tmp, &mApp->worker.model->wp.totalDesorbedMolecules,
                           "%.4g");
        }

        ImGui::EndDisabled();
        //***********************

        {
            const std::string btnText = "Recalc. outgassing";
            ImGui::PlaceAtRegionRight(btnText.c_str(), false);
        }

        if (ImGui::Button("Recalc. outgassing")) // Edit bools storing our
            // window open/close state
            recalcOutg = true;
        simChanged = ImGui::Checkbox(
                "Enable low flux mode",
                &lowFluxMode); // Edit bools storing our window open/close state
        ImGui::SameLine();
        ImGui::HelpMarker(
                "Low flux mode helps to gain more statistics on low pressure "
                "parts of the system, at the expense\n"
                "of higher pressure parts. If a traced particle reflects from a "
                "high sticking factor surface, regardless of that probability,\n"
                "a reflected test particle representing a reduced flux will "
                "still be traced. Therefore test particles can reach low flux "
                "areas more easily, but\n"
                "at the same time tracing a test particle takes longer. The "
                "cutoff ratio defines what ratio of the originally generated "
                "flux\n"
                "can be neglected. If, for example, it is 0.001, then, when "
                "after subsequent reflections the test particle carries less "
                "than 0.1%\n"
                "of the original flux, it will be eliminated. A good advice is "
                "that if you'd like to see pressure across N orders of "
                "magnitude, set it to 1E-N");
        if (!lowFluxMode) {
            ImGui::BeginDisabled();
        }
        simChanged = ImGui::InputRightSide(
                "Cutoff ratio", &lowFluxCutoff,
                "%.2e"); // Edit bools storing our window open/close state
        if (!lowFluxMode) {
            ImGui::EndDisabled();
        }

        {
            ImGui::PlaceAtRegionCenter("Apply above settings");
            if (ImGui::Button("Apply above settings")) {
                simChanged = false;
                mApp->worker.model->wp.gasMass = gasMass;
                mApp->worker.model->wp.enableDecay = enableDecay;
                mApp->worker.model->wp.halfLife = halfLife;
                mApp->worker.model->otfParams.lowFluxMode = lowFluxMode;
                mApp->worker.model->otfParams.lowFluxCutoff = lowFluxCutoff;
            }
        }

#else
//SYNRAD
#endif
        ImGui::PopItemWidth();
        // PopStyleCompact();
        ImGui::EndTable();
    }

    ImGui::Text("Number of CPU cores:     %zd", mApp->numCPU);
    ImGui::Text("Number of subprocesses:  ");
    ImGui::SameLine();
    /*bool nbProcChanged = false;
    int nbProc = mApp->worker.GetProcNumber();*/
    ImGui::SetNextItemWidth(50.0f);
    ImGui::DragInt("", &nbProc, 1, 0, MAX_PROCESS, "%d",
                   ImGuiSliderFlags_AlwaysClamp);

    ImGui::SameLine();
    if (ImGui::Button("Apply and restart processes"))
        nbProcChanged = true;
    {
        ImGui::PlaceAtRegionRight("Change MAX desorbed molecules", true);
        if (ImGui::Button("Change MAX desorbed molecules"))
            ImGui::OpenPopup("Edit MAX");
    }
    // Always center this window when appearing
    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing,
                            ImVec2(0.5f, 0.5f));

    if (ImGui::BeginPopupModal("Edit MAX", NULL,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
        // static bool initMax = false;
        static double maxDes = mApp->worker.model->otfParams.desorptionLimit;
        /*if(!initMax) {
            maxDes = mApp->worker.model->otfParams.desorptionLimit; // use
        double function to allow exponential format initMax = true;
        }*/
        ImGui::InputDouble("Desorption max (0=>endless)", &maxDes, 0.0f, 0.0f,
                           "%.0f");
        ImGui::Separator();

        // static int unused_i = 0;
        // ImGui::Combo("Combo", &unused_i, "Delete\0Delete harder\0");

        /*static bool dont_ask_me_next_time = false;
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, 0));
        ImGui::Checkbox("Don't ask me next time", &dont_ask_me_next_time);
        ImGui::PopStyleVar();*/

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            mApp->worker.model->otfParams.desorptionLimit = maxDes;
            // initMax = false;
            changeDesLimit = true;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SetItemDefaultFocus();
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

}