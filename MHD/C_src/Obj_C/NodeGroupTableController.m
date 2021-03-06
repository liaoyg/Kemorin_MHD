//
//  NodeGroupTableController.m
//  Kemoview_Cocoa
//
//  Created by Hiroaki Matsui on 10/09/29.
//  Copyright 2010 Department of Geophysical Sciences, University of Chicago. All rights reserved.
//

#import "NodeGroupTableController.h"
#include "kemoviewer.h"


@implementation NodeGroupTableController
- (id) init
{
	NumNodeGroup = 0;
	NodeGroupDisplayNames= [[NSMutableArray alloc] init];
	NodeGroupDisplayNodeFlags=  [[NSMutableArray alloc] init];
	
	return self;
}

- (id) dealloc
{
	[NodeGroupDisplayNames release];
	[NodeGroupDisplayNodeFlags  release];
	
	[super dealloc];
	return self;
}

- (IBAction) ShowAllNodeGroupAction:(id)pId{
	int i;
	
	NSLog(@"selectedNodeGroupObjectType %@", selectedNodeGroupObjectType);
	if([selectedNodeGroupObjectType isEqualToString:@"NodGrpNode"]) {
		NSLog(@"Set all NodGrpNode");
		[NodeGroupDisplayNodeFlags removeAllObjects];
		for(i=0;i<NumNodeGroup;i++){
			[NodeGroupDisplayNodeFlags addObject:[[NSNumber alloc ] initWithInt:1] ];
			NSLog(@"Set all nodes %d", i);
			kemoview_set_draw_nodgrp_node(IONE,i);
		}
	}

    [self UpdateNodeTable];
	[_kemoviewer UpdateImage];
}

- (IBAction) HideAllNodeGroupAction:(id)pId
{
	int i;
	
	if([selectedNodeGroupObjectType isEqualToString:@"NodGrpNode"]) {
		[NodeGroupDisplayNodeFlags removeAllObjects];
		for(i=0;i<NumNodeGroup;i++){
			[NodeGroupDisplayNodeFlags addObject:[[NSNumber alloc ] initWithInt:0] ];
			kemoview_set_draw_nodgrp_node(IZERO,i);
		}
	}	

    [self UpdateNodeTable];
	[_kemoviewer UpdateImage];
}

- (int)numberOfRowsInTableView:(NSTableView *)aTableView
{
    return NumNodeGroup;
}

- (id)tableView:(NSTableView *)aTableView
objectValueForTableColumn:(NSTableColumn *)aTableColumn
			row:(int)rowIndex
{
	if([[aTableColumn identifier] isEqualToString:@"NodGrpNode"]) {
        return [NodeGroupDisplayNodeFlags objectAtIndex:rowIndex];
    }
    else if([[aTableColumn identifier] isEqualToString:@"NodGrpName"]) {
		return [NodeGroupDisplayNames objectAtIndex:rowIndex];
	}
    return nil;
}

- (void)tableView:(NSTableView *)aTableView
   setObjectValue:(id)object 
   forTableColumn:(NSTableColumn *)tableColumn 
			  row:(int)rowIndex;
{
    id	identifier;
	
    identifier = [tableColumn identifier];
    if([identifier isEqualToString:@"NodGrpNode"]) {
		kemoview_set_draw_nodgrp_node([object intValue],rowIndex);
		[NodeGroupDisplayNodeFlags replaceObjectAtIndex:rowIndex withObject:object];
	}

	[_kemoviewer UpdateImage];
}

- (void)tableView:(NSTableView *)aTableView didClickTableColumn:(NSTableColumn *)tableColumn
{
	selectedNodeGroupObjectType = [tableColumn identifier];
	return;
}

- (void) UpdateNodeTable
{
	int i, iflag;
	char name[4096];
	NSString *stname;
	
//	printf("Update Node group map\n");
	[NodeGroupDisplayNames removeAllObjects];
	[NodeGroupDisplayNodeFlags removeAllObjects];
	NumNodeGroup = kemoview_get_num_node_grp();
	for(i=0;i<NumNodeGroup;i++){
		kemoview_get_node_grp_name(name,i);
		stname = [[NSString alloc] initWithUTF8String:name];
		[NodeGroupDisplayNames      addObject:stname];
		iflag = kemoview_get_draw_nodgrp_node(i);
		[NodeGroupDisplayNodeFlags  addObject:[[NSNumber alloc ] initWithInt:iflag] ];
		[stname release];
	}
	[_nodeTableView reloadData];
	
}

- (IBAction)ChooseNodeGrpNodeColorAction:(id)sender;
{
	NSInteger tag = [[_NodeGrpNodeColorItem selectedCell] tag];
	kemoview_set_node_grp_color_flag(tag);

	[_kemoviewer UpdateImage];
}
- (IBAction)SetNodeGrpNodeColorAction:(id)sender
{
	CGFloat redBG, greenBG, blueBG, opacityBG;
	GLfloat colorcode4[4];
	nsNodeGrpNodeColor = [nodeGrpNodeColorWell color];
	[nsNodeGrpNodeColor getRed:&redBG green:&greenBG blue:&blueBG alpha:&opacityBG ];
	colorcode4[0] =  (GLfloat) redBG;
	colorcode4[1] =  (GLfloat) greenBG;
	colorcode4[2] =  (GLfloat) blueBG;
	colorcode4[3] =  (GLfloat) opacityBG;
	kemoview_set_node_grp_color_code(colorcode4);
	
	[_kemoviewer UpdateImage];
}
@end
