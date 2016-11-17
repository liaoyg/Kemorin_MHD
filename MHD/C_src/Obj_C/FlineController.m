//
//  FlineController.m
//  Kemoview_Cocoa
//
//  Created by Hiroaki Matsui on 11/08/17.
//  Copyright 2011 Dept. of Earth and Planetary Science, UC Berkeley. All rights reserved.
//

#import "FlineController.h"
#include "kemoviewer.h"


@implementation FlineController

@synthesize FlineOpenFilename;
@synthesize FlineWindowlabel;

@synthesize DrawFlineFlag;
@synthesize FlineMinimumValue;
@synthesize FlineMaximumValue;
@synthesize FlineDisplayMinimum;
@synthesize FlineDisplayMaximum;
@synthesize Flinetype;
@synthesize FlineThickFactor;
@synthesize FlineThickDigit;

- (id)init;
{
	FlineThickFactor = 1;
	FlineThickDigit = -2;
	self.FlineWindowlabel = [NSString stringWithFormat:@"Fieldline View"];
	
	FlineNumberOfComponent =[[NSMutableArray alloc] init];
	FlineFieldName =        [[NSMutableArray alloc] init];
	FlineMinimum =          [[NSMutableArray alloc] init];
	FlineMaximum =          [[NSMutableArray alloc] init];
	
	FieldlineFlag =  [NSNumber alloc];
	
	FlineDrawFieldId =     [NSNumber alloc];
	FlineDrawComponentId = [NSNumber alloc];
	
	FieldlineColor =      [NSNumber alloc];
	
	FlineMinimumRange = [NSNumber alloc];
	FlineMaximumRange = [NSNumber alloc];
	
	return self;
}

- (id)dealloc
{
	
	[FlineNumberOfComponent dealloc];
	[FlineFieldName         dealloc];
	[FlineMinimum           dealloc];
	[FlineMaximum           dealloc];
	
	[FieldlineFlag    dealloc];
	
	[FlineDrawFieldId      dealloc];
	[FlineDrawComponentId  dealloc];
	
	[FieldlineColor      dealloc];

	[FlineMinimumRange dealloc];
	[FlineMaximumRange dealloc];
	
    [super dealloc];
	return self;
}

- (id) CopyFlineDisplayFlagsFromC
{
	int i, iflag;
	double minmax;
	char name[4096];
	NSString *stname;
	NSNumber *stnum;
	
	FlineNumberOfField =  send_nfield_fline();
	FlineTotalComponent = send_ncomptot_fline();
	
	FlineThickFactor = 1;
	FlineThickDigit = -2;
	set_to_fline_thickness(0.01);


	[FlineFieldName removeAllObjects];	
	[FlineNumberOfComponent removeAllObjects];
	[FlineMinimum removeAllObjects];
	[FlineMaximum removeAllObjects];
	for(i = 0; i < FlineNumberOfField; i++){
		send_fline_data_name(name,i);
		stname = [[NSString alloc] initWithUTF8String:name];
		[FlineFieldName      addObject:stname];
		[stname release];
		
		iflag = send_ncomp_fline(i);
		stnum = [[NSNumber alloc] initWithInt:iflag];
		[FlineNumberOfComponent addObject:stnum];
		[stnum release];	
	}
	for(i = 0; i < FlineTotalComponent; i++){
		minmax = send_fline_data_min(i);
		stnum = [[NSNumber alloc] initWithDouble:minmax];
		[FlineMinimum      addObject:stnum];
		[stnum release];	
		
		minmax = send_fline_data_max(i);
		stnum = [[NSNumber alloc] initWithDouble:minmax];
		[FlineMaximum      addObject:stnum];
		[stnum release];	
	}

	for(i = 0; i < FlineNumberOfField; i++){
		NSLog(@"FlineNumberOfComponent = %d %@\n",i, [FlineNumberOfComponent objectAtIndex: i]);
		NSLog(@"FlineFieldName = %d %@\n",i, [FlineFieldName objectAtIndex: i]);
	}
	/*
	NSLog(@"FlineTotalComponent = %d\n",FlineTotalComponent);
	for(i = 0; i < FlineTotalComponent; i++){
		NSLog(@"FlineMinimum = %d %@\n",i, [FlineMinimum objectAtIndex: i]);
		NSLog(@"FlineMaximum = %d %@\n",i, [FlineMaximum objectAtIndex: i]);
	}
	*/

	[self SetFlineFieldMenu];
	[self SetFlineComponentMenu:0];
	return self;
}

- (void) OpenFieldlineFile:(NSString*) fieldlineFilehead{
	int id_viewtype;
	
	self.FlineWindowlabel = [NSString stringWithFormat:@"Fieldline:%@",
							 [[fieldlineFilehead lastPathComponent] stringByDeletingPathExtension]];

	self.DrawFlineFlag = send_iflag_draw_fline();
	[self CopyFlineDisplayFlagsFromC];
	//		self.EvolutionStartStep = [[FlineOpenFilehead pathExtension] intValue];
	//		self.EvolutionEndStep =    self.EvolutionStartStep;
	
	id_viewtype = send_iflag_view_type();
	[_kemoviewControl SetViewTypeMenu:id_viewtype];
	
	
	[_kemoviewer setViewerType:VIEW_3D];
	[_kemoviewer UpdateImage];
};

- (IBAction) DrawFlineFile:(id)pId{
	NSArray *flineFileTypes = [NSArray arrayWithObjects:@"inp",@"vtk",@"gz",@"INP",@"VTK",@"GZ",nil];
	NSOpenPanel *flineOpenPanelObj	= [NSOpenPanel openPanel];
	[flineOpenPanelObj setTitle:@"Choose field line data"];
    [flineOpenPanelObj setAllowedFileTypes:flineFileTypes];
	NSInteger FlineOpenInteger	= [flineOpenPanelObj runModal];
	
	if(FlineOpenInteger == NSFileHandlingPanelOKButton){
		FlineOpenDirectory = [[flineOpenPanelObj directoryURL] path];
		self.FlineOpenFilename =  [[flineOpenPanelObj URL] path];
		FlineOpenFileext =   [self.FlineOpenFilename pathExtension];
		FlineOpenFilehead =  [self.FlineOpenFilename stringByDeletingPathExtension];
		// NSLog(@"PSF file directory = %@",FlineOpenDirectory);
		// NSLog(@"PSF file name =      %@",FlineOpenFilename);
		// NSLog(@"PSF file header =    %@",FlineOpenFilehead);
		// NSLog(@"self.FlineWindowlabel = %@",self.FlineWindowlabel);
		
		if([FlineOpenFileext isEqualToString:@"gz"] || [FlineOpenFileext isEqualToString:@"GZ"]){
			FlineOpenFileext =    [FlineOpenFilehead pathExtension];
			FlineOpenFilehead =   [FlineOpenFilehead stringByDeletingPathExtension];
		};
	
		int iflag_datatype =  kemoview_open_data_glut([self.FlineOpenFilename UTF8String]);
		if(iflag_datatype == IFLAG_LINES) [self OpenFieldlineFile:(NSString *)FlineOpenFilehead];
	};	

}

- (IBAction) CloseFlineFile:(id)pId{

	close_fline_view();
	self.DrawFlineFlag = send_iflag_draw_fline();
	[self CopyFlineDisplayFlagsFromC];
	
	[_kemoviewer UpdateImage];
	
}

- (IBAction) FlineFieldAction:(id)sender
{	
	int isel, iplotted;
	isel = [_FlineFieldMenu indexOfSelectedItem];
	[self SetFlineComponentMenu:isel];
	
	set_fline_color_field_flag(isel);
	
	
	iplotted = send_icomp_draw_fline();
	self.FlineMinimumValue =   send_fline_data_min(iplotted);
	self.FlineMaximumValue =   send_fline_data_max(iplotted);
	self.FlineDisplayMinimum = [[FlineMinimum objectAtIndex:iplotted] doubleValue];
	self.FlineDisplayMaximum = [[FlineMaximum objectAtIndex:iplotted] doubleValue];

	[_kemoviewer UpdateImage];
}

- (IBAction) FlineComponentAction:(id)sender
{	
	int iplotted;
	
	set_fline_color_comp_flag((int) [_FlineComponentMenu indexOfSelectedItem]);
	
	iplotted = send_icomp_draw_fline();
	self.FlineMinimumValue =   send_fline_data_min(iplotted);
	self.FlineMaximumValue =   send_fline_data_max(iplotted);
	self.FlineDisplayMinimum = [[FlineMinimum objectAtIndex:iplotted] doubleValue];
	self.FlineDisplayMaximum = [[FlineMaximum objectAtIndex:iplotted] doubleValue];	

	[_kemoviewer UpdateImage];
}

- (void) SetFlineFieldMenu{
	int i;
	
	if(FlineNumberOfField > 0) [_FlineFieldMenu removeAllItems];
	if(FlineNumberOfField < 1){
		[_FlineFieldMenu addItemWithTitle:@"No Field"];
	} else {
		for(i = 0; i < FlineNumberOfField; i++){
			[_FlineFieldMenu addItemWithTitle:[FlineFieldName objectAtIndex: i]];
		};
	}
}

- (void) SetFlineComponentMenu:(int)isel{
	int iplotted;

	[_FlineComponentMenu removeAllItems];
	// NSLog ([NSString stringWithFormat:@"component %@\n", [FlineNumberOfComponent objectAtIndex:isel]]);	
	
	if (FlineNumberOfField < 1) {
		[_FlineComponentMenu addItemWithTitle:@"No Field"];
	} else {
		if([[FlineNumberOfComponent objectAtIndex:isel] intValue] == 1){
			[_FlineComponentMenu addItemWithTitle:@"Scalar"];
		}
		else if([[FlineNumberOfComponent objectAtIndex:isel] intValue] == 6){
			[_FlineComponentMenu addItemWithTitle:@"xx"];
			[_FlineComponentMenu addItemWithTitle:@"xy"];
			[_FlineComponentMenu addItemWithTitle:@"xz"];
			[_FlineComponentMenu addItemWithTitle:@"yy"];
			[_FlineComponentMenu addItemWithTitle:@"yz"];
			[_FlineComponentMenu addItemWithTitle:@"zz"];
		}
		else if([[FlineNumberOfComponent objectAtIndex:isel] intValue] == 3){
			NSInteger charalen = [[FlineFieldName objectAtIndex:isel] length];
			if(charalen > 4){
				NSString *stname = [[FlineFieldName objectAtIndex:isel] substringFromIndex:charalen-4];
				// NSLog ([NSString stringWithFormat:@"end is %@\n",stname ]);
				if([stname compare:@"_sph"] == NSOrderedSame){
					[_FlineComponentMenu addItemWithTitle:@"r"];
					[_FlineComponentMenu addItemWithTitle:@"θ"];
					[_FlineComponentMenu addItemWithTitle:@"φ"];
				} else if([stname compare:@"_cyl"] == NSOrderedSame){
					[_FlineComponentMenu addItemWithTitle:@"s"];
					[_FlineComponentMenu addItemWithTitle:@"φ"];
					[_FlineComponentMenu addItemWithTitle:@"z"];
				} else{
					[_FlineComponentMenu addItemWithTitle:@"x"];
					[_FlineComponentMenu addItemWithTitle:@"y"];
					[_FlineComponentMenu addItemWithTitle:@"z"];
				}
			} else {
				[_FlineComponentMenu addItemWithTitle:@"x"];
				[_FlineComponentMenu addItemWithTitle:@"y"];
				[_FlineComponentMenu addItemWithTitle:@"z"];
			}
		}

		iplotted = send_icomp_draw_fline();
		self.FlineMinimumValue =   send_fline_data_min(iplotted);
		self.FlineMaximumValue =   send_fline_data_max(iplotted);
		self.FlineDisplayMinimum = [[FlineMinimum objectAtIndex:iplotted] doubleValue];
		self.FlineDisplayMaximum = [[FlineMaximum objectAtIndex:iplotted] doubleValue];
	}
	return;
}

- (IBAction)ChooseFieldlineColorAction:(id)sender;
{
	NSInteger tag = [[FieldlineColorItem selectedCell] tag];
	set_to_fieldline_color((int) tag);
	
	[_kemoviewer UpdateImage];
}

- (IBAction) ShowFlineRange:(id)pSender {
	input_fline_linear_colormap(self.FlineDisplayMinimum, self.FlineDisplayMaximum);
	[_kemoviewer UpdateImage];
}

- (IBAction)ChooseFieldlineTypeAction:(id)sender;
{
	self.Flinetype = [[_flinetype_matrix selectedCell] tag];
	set_to_fline_type_flag((int) self.Flinetype);
	
	[_kemoviewer UpdateImage];
}

- (IBAction)SetFieldlineThicknessAction:(id)sender;
{
	float fieldlineThickness;
	int IntFlineThickFactor, IntFlineThickDigit;
	int i;
	
	IntFlineThickFactor = (int) FlineThickFactor;
	IntFlineThickDigit = (int) FlineThickDigit;
	
	fieldlineThickness = (double) IntFlineThickFactor;
	if(IntFlineThickDigit < 0){
		for(i=0;i< (-IntFlineThickDigit);i++) fieldlineThickness = fieldlineThickness / 10.0;
	}
	else {
		for(i=0;i<IntFlineThickDigit;i++) fieldlineThickness = fieldlineThickness * 10.0;
	}
	
	set_to_fline_thickness((double) fieldlineThickness);
	
	[_kemoviewer UpdateImage];
}

@end
