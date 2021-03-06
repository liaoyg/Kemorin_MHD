//
//  RGBAMapController.m
//  025-NSTableView
//
//  Created by Hiroaki Matsui on 11/08/23.
//  Copyright 2011 Dept. of Earth and Planetary Science, UC Berkeley. All rights reserved.
//

#import "RGBAMapController.h"
#include "kemoviewer.h"


@implementation RGBAMapController
@synthesize DataMaximum;
@synthesize DataMinimum;

- (void)awakeFromNib {
	self.DataMinimum = ZERO;
	self.DataMaximum = ONE;
}

- (void)updateColormapParameter {
	self.DataMinimum = kemoview_get_PSF_color_table_min();
	self.DataMaximum = kemoview_get_PSF_color_table_max();
}

- (IBAction)SetColorMode:(id)pId;
{
	kemoview_set_PSF_color_mode((int) [ColorModeItem indexOfSelectedItem]);
	[_kemoviewer UpdateImage];
}

- (void) SetColormapMinMax{
	[_colorMapObject InitColorTables];
	[_colorMapObject SetColorTables];
	[_opacityMapObject InitOpacityTables];
	[_opacityMapObject SetOpacityTables];
}


- (IBAction) SaveColormapFile:(id)pId;{
	
	NSSavePanel *ColormapSavePanelObj	= [NSSavePanel savePanel];
    [ColormapSavePanelObj beginSheetModalForWindow:window 
                                    completionHandler:^(NSInteger ColrmapSaveInt){
	if(ColrmapSaveInt == NSFileHandlingPanelOKButton){
		
		NSString * ColormapFilename = [[ ColormapSavePanelObj URL] path];
		NSString * ColormapDirectory = [[ ColormapSavePanelObj directoryURL] path];
		NSString * ColormapFilehead = [ ColormapFilename stringByDeletingPathExtension];
		NSLog(@" ColormapFilename = %@",  ColormapFilename);
		NSLog(@" ColormapDirectory = %@", ColormapDirectory);
		NSLog(@" ColormapFilehead = %@",  ColormapFilehead);
		
		kemoview_write_PSF_colormap_file([ColormapFilename UTF8String]);
	};
                                    }];
}

- (IBAction) LoadColormapFile:(id)pId;{
    NSArray *ColormapFileTypes = [NSArray arrayWithObjects:@"dat",@"DAT",nil];
    
    NSOpenPanel *ColormapOpenPanelObj	= [NSOpenPanel openPanel];
    [ColormapOpenPanelObj setTitle:@"Choose colormap data"];
    [ColormapOpenPanelObj setAllowedFileTypes:ColormapFileTypes];
    [ColormapOpenPanelObj beginSheetModalForWindow:window 
                                   completionHandler:^(NSInteger ColormapOpenInteger){
                                       if(ColormapOpenInteger == NSFileHandlingPanelOKButton){
                                           
                                           NSString * ColormapFilename = [[ ColormapOpenPanelObj URL] path];
                                           NSString * ColormapDirectory = [[ ColormapOpenPanelObj directoryURL] path];
                                           NSString * ColormapFilehead = [ ColormapFilename stringByDeletingPathExtension];
                                           NSLog(@" ColormapFilename = %@",  ColormapFilename);
                                           NSLog(@" ColormapDirectory = %@", ColormapDirectory);
                                           NSLog(@" ColormapFilehead = %@",  ColormapFilehead);
                                           
                                           kemoview_read_PSF_colormap_file((char *) [ColormapFilename UTF8String]);
                                           [_kemoviewer UpdateImage];
                                           [_colorMapObject SetColorTables];
                                           [_opacityMapObject SetOpacityTables];
                                           [ColorModeItem selectItemAtIndex:kemoview_get_PSF_color_mode()];
                                       };
                                   }];
    
}

@end
