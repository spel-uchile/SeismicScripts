#!/usr/bin/Python
# -*- coding: utf-8 -*-

# from dbus import exceptions
from datetime import datetime
import dateutil.parser
import datetime
import os

__author__ = 'toopazo'


class ParserUtilities():
    def __init__(self):
        pass

    @staticmethod
    def get_xxx_files(folderpath, extension):
        xxx_files = []
        for file_i in os.listdir(folderpath):
            if os.path.isfile(os.path.join(folderpath, file_i)):
                # filename, file_extension = os.path.splitext(file_i)
                # if file_extension == extension:
                #     xxx_files.append(file_i)
                if extension in file_i:
                    xxx_files.append(file_i)
        xxx_files.sort()
        nfiles = len(xxx_files)
        print "[get_xxx_files] %s \"%s\" files were found at %s" % (nfiles, extension, folderpath)

        return xxx_files

    @staticmethod
    def get_all_files(folderpath):
        xxx_files = []
        for file_i in os.listdir(folderpath):
            if os.path.isfile(os.path.join(folderpath, file_i)):
                # filename, file_extension = os.path.splitext(file_i)
                # if file_extension == extension:
                #     xxx_files.append(file_i)
                xxx_files.append(file_i)
        xxx_files.sort()
        nfiles = len(xxx_files)
        print "[get_xxx_files] %s files were found at %s" % (nfiles, folderpath)

        return xxx_files