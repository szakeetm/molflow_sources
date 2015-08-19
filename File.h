/*
  File:        File.h
  Description: File management class
  Program:     MolFlow
  Author:      R. KERSEVAN / J-L PONS / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef FILERWH
#define FILERWH

#include <stdio.h>
#include <string>
#include "Types.h"

#define READ_BUFFSIZE 4096

// -------------------------------------------------------

class Error {

public:
  Error(const char *message);
  const char *GetMsg();

private:
  char msg[1024];

};

// -------------------------------------------------------

class FileUtils {

public:
  // Utils functions
  static int Exist(const char *fileName);
  static int Exist(std::string fileName);
  static char *GetPath(char *fileName);
  static std::string GetPath(const std::string &str); //Extracts string up to to last "\"
  static std::string GetFilename(const std::string &str); //Extracts string after the last "\"
  static std::string GetExtension(const std::string &str); //Extracts string after the last "."

};

// -------------------------------------------------------

class FileReader {

public:
  // Constructor/Destructor
  FileReader(char *fileName);
  FileReader(std::string fileName);
  ~FileReader();

  char *GetName();

  // Read function
  int IsEof();
  char *ReadLine();
  char *ReadString();
  llong ReadLLong();
  int ReadInt();
  double ReadDouble();
  void ReadKeyword(char *keyword);
  char *ReadWord();
  void JumpSection(char *end);
  void SeekStart();
  bool SeekFor(char *keyword);


  Error MakeError(char *msg);
  int GetCurrentLine();

private:

  void RefillBuffer();
  void ReadChar();
  void JumpSpace();
  FILE *file;
  int curLine;
  char fileName[1024];
  char readBuffer[READ_BUFFSIZE];
  int  nbLeft;
  int  buffPos;
  int  isEof;
  char CurrentChar;

};

// -------------------------------------------------------

class FileWriter {

public:
  // Constructor/Destructor
  FileWriter(char *fileName);
  ~FileWriter();

  char *GetName();

  // Write function
  void WriteLLong(llong v,char *sep=NULL);
  void WriteInt(int v,char *sep=NULL);
  void WriteDouble(double v,char *sep=NULL);
  void Write(const char *s);

private:

  FILE *file;
  char fileName[1024];

};

#endif /* FILERWH */

