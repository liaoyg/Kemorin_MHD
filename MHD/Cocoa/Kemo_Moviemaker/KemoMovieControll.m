//
//  KemoMovieControll.m
//  Kemo_Moviemaker
//
//  Created by Hiroaki Matsui on 10/10/26.
//  Copyright 2010 Department of Geophysical Sciences, University of Chicago. All rights reserved.
//

#import "KemoMovieControll.h"

@implementation KemoMovieControll
@synthesize evolutionCurrentStep;
@synthesize evolutionStartStep;
@synthesize evolutionEndStep;
@synthesize evolutionIncrement;
@synthesize evolutionFPS;

- (id)init;
{
	self.evolutionCurrentStep = 1;
	self.evolutionStartStep = 1;
	self.evolutionEndStep =   1;
	self.evolutionIncrement = 1;
	self.evolutionFPS = 12;
	return self;
}	

-(CGImageRef) NSImageToCGimageref:(NSImage*)nsImage{
//    NSString *path = /* Take a path to load the file ... */
//    NSImage *nsImage = [[NSImage alloc] initWithContentsOfFile:path];
    CGImageRef cgImage;
    if (nsImage != nil) {
        NSSize size = [nsImage size];
        uint32_t width = (uint32_t) size.width, height = (uint32_t) size.height, components = 4;
        uint8_t *pixels = (uint8_t *) malloc(size.width * size.height * components);
        if (pixels) {
            CGColorSpaceRef colorSpace = CGColorSpaceCreateDeviceRGB();
            CGContextRef bitmapContext = CGBitmapContextCreate(pixels, width, height, 8, components * width, colorSpace, kCGImageAlphaPremultipliedLast);
            NSRect rect = NSMakeRect(0, 0, width, height);
            NSGraphicsContext *graphicsContext = (NSGraphicsContext *) [[NSGraphicsContext currentContext] graphicsPort];
            cgImage = [nsImage CGImageForProposedRect:&rect context:graphicsContext hints:nil];
            CGContextDrawImage(bitmapContext, NSRectToCGRect(rect), cgImage);
            CGContextRelease(bitmapContext);
            CGColorSpaceRelease(colorSpace);
        /* Handle pixels ... */
            free(pixels);
        }
    }
    return cgImage;
}

- (CVPixelBufferRef)pixelBufferFromCGImage:(CGImageRef)cgImage
{
    NSDictionary *options = @{ (NSString *)kCVPixelBufferCGImageCompatibilityKey: @YES,
                               (NSString *)kCVPixelBufferCGBitmapContextCompatibilityKey: @YES, };
    
    CVPixelBufferRef pxbuffer = NULL;
    
    CGFloat width  = CGImageGetWidth(cgImage);
    CGFloat height = CGImageGetHeight(cgImage);
    CVPixelBufferCreate(kCFAllocatorDefault,
                        width,
                        height,
                        kCVPixelFormatType_32ARGB,
                        (__bridge CFDictionaryRef)options,
                        &pxbuffer);
    
    CVPixelBufferLockBaseAddress(pxbuffer, 0);
    void *pxdata = CVPixelBufferGetBaseAddress(pxbuffer);
    
    size_t bitsPerComponent       = 8;
    size_t bytesPerRow            = 4 * width;
    CGColorSpaceRef rgbColorSpace = CGColorSpaceCreateDeviceRGB();
    CGContextRef context = CGBitmapContextCreate(pxdata,
                                                 width,
                                                 height,
                                                 bitsPerComponent,
                                                 bytesPerRow,
                                                 rgbColorSpace,
                                                 (CGBitmapInfo)kCGImageAlphaNoneSkipFirst);
    
    CGContextDrawImage(context, CGRectMake(0, 0, width, height), cgImage);
    CGColorSpaceRelease(rgbColorSpace);
    CGContextRelease(context);
    
    CVPixelBufferUnlockBaseAddress(pxbuffer, 0);
    
    return pxbuffer;
}



-(IBAction) SaveImageEvolution:(id)pSender
{
    NSImage *anImage;
	//	NSLog(@"SaveRotation received message = %@",(NSString*)[pNotification object]);
	NSError *overWriteflag = [[NSError alloc] init];
	
	NSOpenPanel *PsfOpenPanelObj	= [NSOpenPanel openPanel];
	[PsfOpenPanelObj setTitle:@"Choose one of image files"];
	if([PsfOpenPanelObj runModal] == NSModalResponseOK){
/*
        NSURL *PsfURL = [[NSURL alloc] init];
        PsfURL = [PsfOpenPanelObj URL];
        NSLog(@"absoluteString : %@", [PsfURL absoluteString]);
        NSLog(@"absoluteURL : %@", [PsfURL absoluteURL]);
        NSLog(@"baseURL : %@", [PsfURL baseURL]);
        NSLog(@"fragment : %@", [PsfURL fragment]);
        NSLog(@"host : %@", [PsfURL host]);
        NSLog(@"lastPathComponent : %@", [PsfURL lastPathComponent]);
        NSLog(@"parameterString : %@", [PsfURL parameterString]);
        NSLog(@"password : %@", [PsfURL password]);
        NSLog(@"path : %@", [PsfURL path]);
        NSLog(@"pathComponents : %@", [PsfURL pathComponents]);
        NSLog(@"pathExtension : %@", [PsfURL pathExtension]);
        NSLog(@"port : %@", [PsfURL port]);
        NSLog(@"query : %@", [PsfURL query]);
        NSLog(@"relativePath : %@", [PsfURL relativePath]);
        NSLog(@"relativeString : %@", [PsfURL relativeString]);
        NSLog(@"resourceSpecifier : %@", [PsfURL resourceSpecifier]);
        NSLog(@"scheme : %@", [PsfURL scheme]);
        NSLog(@"standardizedURL : %@", [PsfURL standardizedURL]);
        NSLog(@"user : %@", [PsfURL user]);
*/
		imageFileName =  [[PsfOpenPanelObj URL] path];
		imageFileExt =   [imageFileName pathExtension];
		imageFileHead =  [imageFileName stringByDeletingPathExtension];
		imageFileHeadExStep =  [imageFileHead stringByDeletingPathExtension];
//        NSLog(@"PSF file name =      %@",imageFileName);
//		NSLog(@"PSF file header =    %@",imageFileHead);
	}
	else { return;};
	
	NSSavePanel *evolutionImageSavePanelObj	= [NSSavePanel savePanel];
	[evolutionImageSavePanelObj setCanSelectHiddenExtension:YES];	

    NSURL *movieFileURL = [movieFileURL init];
    if([evolutionImageSavePanelObj runModal] == NSModalResponseOK){
        movieFileURL = [evolutionImageSavePanelObj URL];
		movieFileName = [movieFileURL path];
		movieFileHead = [movieFileName stringByDeletingPathExtension];
		movieFileName = [movieFileHead stringByAppendingPathExtension:@"mov"];
		NSError *overWriteflag = [[NSError alloc] init];
	}
	else { return;};

    NSError *error = [NSError init];
    if(overWriteflag!= NULL ){
        NSFileManager *fman = [NSFileManager defaultManager];
        [fman removeItemAtURL:movieFileURL error:&error];
    }

    // Initialize AVKit
    anImage = [[NSImage alloc] initWithContentsOfFile:imageFileName];
    NSLog(@"PSF width =      %f",[anImage size].width);
    NSLog(@"PSF width =      %f",[anImage size].height);
    
    NSDictionary *videosettings = [NSDictionary dictionaryWithObjectsAndKeys:
                                  AVVideoCodecH264, AVVideoCodecKey,
                                  [anImage size].width, AVVideoWidthKey,
                                  [anImage size].height, AVVideoHeightKey,
                                  nil];
    NSDictionary *pixelBufferAttributes = [NSDictionary dictionaryWithObjectsAndKeys:
                                                [NSNumber numberWithInt:kCVPixelFormatType_32ARGB],
                                                 kCVPixelBufferPixelFormatTypeKey, nil];

    [anImage release];
   
    AVAssetWriter *videoWriter = [[AVAssetWriter alloc] initWithURL:movieFileURL fileType:AVFileTypeMPEG4 error:&error];
    AVAssetWriterInput *writerInput
        = [AVAssetWriterInput assetWriterInputWithMediaType:AVMediaTypeVideo outputSettings:videosettings];
    AVAssetWriterInputPixelBufferAdaptor *adaptor
        = [AVAssetWriterInputPixelBufferAdaptor assetWriterInputPixelBufferAdaptorWithAssetWriterInput:writerInput sourcePixelBufferAttributes:pixelBufferAttributes];
    
    writerInput.expectsMediaDataInRealTime = YES;
    [videoWriter addInput:writerInput]; 
    
    // Create a QTMovie with a writable data reference
    [videoWriter startWriting];
    [videoWriter startSessionAtSourceTime:kCMTimeZero];
    // pixel bufferを宣言
    CVPixelBufferRef buffer = NULL;
    
    // 現在のフレームカウント
    int frameCount = 0;
    
    // 各画像の表示する時間
    int durationForEachImage = 1;
    
    [progreessBar setIndeterminate:NO];
    [progreessBar startAnimation:(id)pSender];
	for (self.evolutionCurrentStep = self.evolutionStartStep; 
		 self.evolutionCurrentStep<(self.evolutionEndStep+1);
		 self.evolutionCurrentStep++) {
        if( ((self.evolutionCurrentStep-self.evolutionStartStep)%self.evolutionIncrement) == 0) {
        imageFileHead =  [imageFileHeadExStep stringByAppendingPathExtension:
                        [NSString stringWithFormat:@"%d",(int) self.evolutionCurrentStep]];
        imageFileName =  [imageFileHead stringByAppendingPathExtension:imageFileExt];
        [progreessBar incrementBy:(double)self.evolutionIncrement];
        [progreessBar displayIfNeeded];
        @autoreleasepool {
            if (!adaptor.assetWriterInput.readyForMoreMediaData) {break;}
            CMTime frameTime = CMTimeMake((int64_t)frameCount * self.evolutionFPS * durationForEachImage, self.evolutionFPS);

            NSImage *anImage = [[NSImage alloc] initWithContentsOfFile:imageFileName];
            CGImageRef cgImage = [self NSImageToCGimageref:anImage];
            buffer = [self pixelBufferFromCGImage:cgImage];
            if (![adaptor appendPixelBuffer:buffer withPresentationTime:frameTime]) {
                // Error!
            }
            if (buffer) {CVBufferRelease(buffer);}
            frameCount++;
        }
		};
	};
	[progreessBar setDoubleValue:(double)self.evolutionStartStep];
	[progreessBar displayIfNeeded];
    
    [writerInput markAsFinished];
    [videoWriter endSessionAtSourceTime:CMTimeMake((int64_t)(frameCount - 1) * self.evolutionFPS * durationForEachImage, self.evolutionFPS)];
    
    [videoWriter finishWritingWithCompletionHandler:^{
        // Finish!
    }];
    
    CVPixelBufferPoolRelease(adaptor.pixelBufferPool);
	return;
}

- (IBAction)SetEvolutionSteps:(id)pSender{
	NSLog(@"start: %d", (int) self.evolutionStartStep);
	NSLog(@"end: %d", (int) self.evolutionEndStep);
	NSLog(@"increment: %d", (int) self.evolutionIncrement);
	NSLog(@"FPS: %d", (int) self.evolutionFPS);
}

@end
